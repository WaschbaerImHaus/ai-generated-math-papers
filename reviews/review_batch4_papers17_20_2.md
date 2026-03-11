# Gutachter-Bericht: Batch 4 — Papers 17–20
## Review-Datum: 2026-03-11  |  Build 6

---

## Überblick

Batch 4 umfasst vier Papers zur **analytischen Zahlentheorie** (Kreismethode-Serie):

| Paper | Titel (EN) | EN | DE |
|-------|-----------|----|----|
| 17 | Hardy–Littlewood Circle Method | paper17_circle_method.tex | paper17_kreismethode_de.tex |
| 18 | Vinogradov's Three-Prime Theorem | paper18_vinogradov_three_primes.tex | paper18_vinogradov_drei_primzahlen_de.tex |
| 19 | Waring's Problem | paper19_waring_problem.tex | paper19_waringsches_problem_de.tex |
| 20 | Goldbach Singular Series | paper20_goldbach_singular_series.tex | paper20_goldbach_singulaere_reihe_de.tex |

---

## Paper 17 — Kreismethode (EN + DE)

### Formale Korrektheit
**EN:** LaTeX-Struktur korrekt. Alle Theorem-Environments definiert und korrekt verwendet.
**Hinweis:** Drei `\DeclareMathOperator`-Deklarationen (`\ord`, `\li`, `\sgn`) werden im Paper nie verwendet. LaTeX-Warnung, kein Fehler.

**DE:** Korrekte Präambel mit `inputenc`, `ngerman`, deutschen Theorem-Names. Vollständig kompilierbar.

### Mathematische Schlüssigkeit
- Weyl-Summen (Def. 2.1), Gauß-Summen (Def. 2.3), Farey-Zerlegung (Def. 3.1–3.2): korrekt.
- Lemma 2.2 (geometrische Reihe): Beweis korrekt, $|\sin(\pi x)| \geq 2\|x\|$ korrekt verwendet.
- Haupt-/Nebenbogen-Definition (Def. 4.1) und Voraussetzung $Q^2 \leq N$: korrekt.
- Lemma 5.1 (Hauptbogen-Approximation): Fehlerterm im Beweis zunächst als $O(q)$ hergeleitet, dann auf $O(q^{1+\varepsilon})$ verbessert via Verweis auf Vaughan — korrekt.
- Satz 5.6 (Hauptbogen-Beitrag): Beweisskizze korrekt, Fehlerterme konsistent.
- Weylsche Ungleichung (Satz 6.1) und Korollar 6.2 (Nebenbogen): korrekt.
- Hauptsatz 7.1 (Hardy–Littlewood–Vinogradov): Beweisstruktur korrekt (Parseval + Nebenbogenabschätzung).

### Vollständigkeit
- EN enthält alle Kernbestandteile: Exponentialsummen, Farey-Zerlegung, Haupt-/Nebenbögen, singuläre Reihe & Integral, Hauptsatz, Goldbach-Anwendung, historische Anmerkungen.
- **EN vs. DE-Abgleich:** DE enthält einen **eigenständigeren Beweis** des Farey-Partition-Lemmas (explizites Disjunktheitsargument via Bogenbreiten) gegenüber dem EN-Verweis auf Hardy–Wright. Das ist eine Verbesserung, kein Mangel.
- DE enthält in der historischen Sektion keinen expliziten Verweis auf das Bourgain–Demeter–Guth-Theorem (EN Abschnitt 8). Inhaltlich kein Verlust — für DE ist das angemessen kompakt.
- Alle Definitionen, Lemmata, Sätze, Beweise oder Beweisszkizzen vorhanden.

### Urteil Paper 17
**EN:** ✅ DRUCKREIF — keine Fehler.
**DE:** ✅ DRUCKREIF — keine Fehler.

---

## Paper 18 — Vinogradov-Drei-Primzahlen-Satz (EN + DE)

### Formale Korrektheit
**EN:** Alle Environments korrekt. `\DeclareMathOperator{\Li}{li}` wird nicht verwendet — LaTeX-Warnung.
**DE:** Korrekte Präambel. Vollständig kompilierbar.

### Mathematische Schlüssigkeit

#### ⚠️ BUG-B4-P18-EN-001 — Falsches Vorzeichen in Remark 4.3 (EN)
**Schwere:** Mittel (Inkonsistenz im Beweis)
**Ort:** EN, Remark 4.3 ("Even $n$ versus odd $n$"), Zeile ~292
**Befund:** Das Remark schreibt für gerades $n$:
$$1 + \frac{c_2(n)}{(2-1)^3} = 1 + c_2(n) = 1 + (2-1) = 2$$
Das Vorzeichen ist **falsch**. Die im Euler-Produkt (Beweis zu Theorem 4.1) korrekt hergeleitete Formel lautet:
$$1 + \left(\frac{\mu(p)}{\phi(p)}\right)^3 c_p(n) = 1 - \frac{c_p(n)}{(p-1)^3}$$
Für $p=2$, gerades $n$: $c_2(n) = 1$, daher Faktor $= 1 - 1 = 0$ ✓ (ergibt $\mathfrak{S}(n) = 0$ für gerades $n$).
Der Faktor $(\mu(2)/\phi(2))^3 = (-1)^3 = -1$, nicht $+1$. Das Remark verwendet fälschlich $+1/(p-1)^3$ statt $-1/(p-1)^3$.
**Folge:** Das Remark kommt auf Faktor $= 2$ statt $= 0$, widerspricht damit der sofort danach folgenden Behauptung $\mathfrak{S}(n) = 0$ und ist intern inkonsistent. Die nachfolgende alternative Erklärung via $\mu(2)^3 e(-n/2)$ ist korrekt, macht das Remark aber noch unübersichtlicher.
**Korrektur:** Das $+$-Zeichen in `$1 + c_2(n)/(2-1)^3$` in `$1 - c_2(n)/(2-1)^3$` ändern.

#### ⚠️ BUG-B4-P18-DE-001 — Gleiches Vorzeichenproblem in Bemerkung (DE)
**Schwere:** Mittel (gleiche Inkonsistenz)
**Ort:** DE, Bemerkung [Gerades $n$ versus ungerades $n$], Zeile ~200
**Befund:** Schreibt `$1 + c_2(n)/(2-1)^3 = 1 + c_2(n)$` statt `$1 - c_2(n)/(2-1)^3$`.
Behauptet danach korrekt $\mathfrak{S}(n) = 0$ für gerades $n$, aber der gezeigten Rechnung widersprechend.
**Korrektur:** Analog zu EN: Vorzeichen korrigieren.

#### ⚠️ BUG-B4-P18-DE-002 — Unvollständiger Beweis des singulären Integrals (DE)
**Schwere:** Gering (inhaltliche Lücke)
**Ort:** DE, Satz 5.1 (Singuläres Integral), Beweis
**Befund:** Der Beweis schreibt:
$$\mathfrak{J}(n) = \text{Fläche}\{t_1+t_2+t_3=n,\; t_i \in [2,N]^3\} = \frac{(n-6)^2}{2} + O(n) \sim \frac{n^2}{2}$$
Die Formel $(n-6)^2/2$ gilt nur für den Fall $n \leq N$ (kleines $n$, Simplex nicht abgeschnitten). Die EN-Version erläutert korrekt drei Fälle: $n \leq N$, $N \leq n \leq 2N$, $2N \leq n \leq 3N$ — der dominante Term ist immer $n^2/2$, aber nur mit Fallunterscheidung korrekt begründet.
**Korrektur:** Fallunterscheidung aus EN ergänzen.

### Vollständigkeit
- EN: sehr vollständig. Euler-Produkt, untere Schranke für $\mathfrak{S}(n)$, Helfgott-Abschnitt, Paritätsproblem — alle Bestandteile vorhanden.
- **EN vs. DE-Abgleich:** DE enthält einen eigenständigen Abschnitt "Das Paritätsproblem" (Kap. 8) der in EN in Remark 7.3 eingebettet ist — strukturell gleichwertig.
- DE fehlt die ausführliche Berechnung des Hauptbogen-Integrals aus EN-Abschnitt 3.2 (Subsection "Computing the major arc integral"). Die Berechnung wird im DE-Beweis von Satz 7.1 stattdessen kompakter durchgeführt — inhaltlich äquivalent.

### Urteil Paper 18
**EN:** ⚠️ ÜBERARBEITUNG EMPFOHLEN — BUG-B4-P18-EN-001 (Vorzeichen in Remark 4.3).
**DE:** ⚠️ ÜBERARBEITUNG EMPFOHLEN — BUG-B4-P18-DE-001 (Vorzeichen), BUG-B4-P18-DE-002 (unvollständiger Beweis Satz 5.1).

---

## Paper 19 — Waringsches Problem (EN + DE)

### Formale Korrektheit
**EN:** Korrekte LaTeX-Struktur. `\newtheorem{proposition}` definiert, nicht verwendet — kein Fehler.
**DE:** Korrekte Präambel. Vollständig.

### Mathematische Schlüssigkeit
- Definitionen g(k), G(k): korrekt.
- Lagranges Vier-Quadrate-Satz: korrekt (Mod-8-Argument für $G(2) \geq 4$, Euler-Identität).
- Hilbert–Waring: korrekte Beweisidee (Hilberts Identität, Dickson-Vereinfachung).
- Hardy–Littlewood-Formel (Satz 4.1 EN / Satz 4.1 DE): korrekt, Verweis auf Paper 17.
- Vinogradov MVT (EN Satz 5.1/5.2): korrekt.
- Wooleys effizientes Kongruenzrechnen (EN Satz 6.1): Der Formalismus des induktiven Schritts ist korrekt dargestellt.
- Formel für $g(k) = 2^k + q - 2$ mit $q = \lfloor(3/2)^k\rfloor$: korrekt.

#### ⚠️ BUG-B4-P19-EN-001 — Falsche Querverweise auf Paper 17 (EN)
**Schwere:** Gering (Referenzfehler)
**Ort:** EN, Theorem 4.1 (Hardy–Littlewood-Formel), Zeile ~206–208
**Befund:** Das Theorem verweist auf "Paper~17, Definition~5.1" und "Paper~17, Lemma~5.3". In Paper 17 EN ist die Nummerierung des Abschnitts 5:
- 5.1 = Lemma (Hauptbogen-Approximation) ← **kein** Definition
- 5.2 = Definition [Singuläre Reihe] ← das gemeinte Objekt
- 5.3 = Remark ← **kein** Lemma
- 5.5 = Lemma [Singuläres Integral, Abschätzung] ← das gemeinte Objekt

Korrekt: "Paper~17, Definition~5.2" und "Paper~17, Lemma~5.5".

#### ⚠️ BUG-B4-P19-DE-001 — Fehlendes Korollar "Weyl-Schranke für G(k)" (DE)
**Schwere:** Mittel (inhaltliche Lücke)
**Ort:** DE, Abschnitt 4 (Kreismethode)
**Befund:** EN enthält Corollary 4.1 mit Beweis: $G(k) \leq (k-2)2^{k-1}+5$ aus Weyls Ungleichung. DE hat nur den Hauptsatz (Hardy–Littlewood-Formel) ohne dieses Korollar.
**Korrektur:** Korollar und Beweis analog EN ergänzen.

#### ⚠️ BUG-B4-P19-DE-002 — Korollar (Wooley-Schranke) ohne Beweis (DE)
**Schwere:** Gering
**Ort:** DE, Abschnitt 5 (Winogradow-Mittelwertsatz), Korollar 5.3
**Befund:** EN enthält Corollary 5.3 (Wooley's $G(k)$ bound) mit einem Beweis (~10 Zeilen). DE hat nur die Behauptung.
**Korrektur:** Beweisszkizze aus EN adaptieren.

#### ⚠️ BUG-B4-P19-DE-003 — Kein formaler Satz für effizientes Kongruenzrechnen (DE)
**Schwere:** Mittel (inhaltliche Lücke)
**Ort:** DE, Abschnitt 6 (Wooleys effizientes Kongruenzrechnen)
**Befund:** EN enthält Definition 6.1 ($p$-adisches Kongruenzrechnen), Theorem 6.1 (Efficient Congruencing, Wooley 2016) mit formaler Aussage und Beweisidee. DE ersetzt dies durch zwei `\bemerkung`-Umgebungen ohne formalen Satz. Kein Theorem-Statement mit präziser Formel für $J_{s,k}$.
**Korrektur:** Theorem 6.1 aus EN adaptieren (DE-Äquivalent mit formaler Formel).

### Vollständigkeit (Tabellen)
- DE-Tabelle für $g(k)$ enthält nur $k = 2, \ldots, 6$ (EN: bis $k = 8$). Geringfügig.
- Beide Tabellen für $G(k)$ (EN und DE) korrekt und konsistent.

### Urteil Paper 19
**EN:** ⚠️ ÜBERARBEITUNG EMPFOHLEN — BUG-B4-P19-EN-001 (Querverweise).
**DE:** ⚠️ ÜBERARBEITUNG EMPFOHLEN — BUG-B4-P19-DE-001, -002, -003 (fehlende Inhalte).

---

## Paper 20 — Goldbach'sche singuläre Reihe (EN + DE)

### Formale Korrektheit
**EN:** Korrekte LaTeX-Struktur. `\DeclareMathOperator{\Li}{li}` nicht verwendet — LaTeX-Warnung.
**DE:** Korrekte Präambel. Vollständig kompilierbar.

### Mathematische Schlüssigkeit

#### Euler-Produkt für $\mathfrak{S}(n)$ (Theorem 2.1 EN / Satz 2.1 DE)
**Formel:** $\mathfrak{S}(n) = 2C_2 \prod_{p \mid n,\, p \geq 3} \frac{p-1}{p-2}$
Beweis durch Faktorisierung der Ramanujan-Summe: korrekt. Die drei Fälle ($p=2$, $p \geq 3$ mit $p \mid n$, $p \geq 3$ mit $p \nmid n$) korrekt berechnet. Konvergenz des Euler-Produkts: korrekt.

**DE-Beweis:** Kompakter, aber korrekt — zeigt $\mathfrak{S}(n) = 2C_2 \prod_{p \geq 3, p \mid n} (p-1)/(p-2)$ durch direktes algebraisches Umformen des Produkts. ✓

#### Positivität (Theorem 3.1 EN / Satz 3.1 DE)
Korrekt: alle Faktoren $> 0$ für gerades $n \geq 4$.

#### ⚠️ BUG-B4-P20-EN-001 — Unvollständiger Satz in Remark (Odd n) (EN)
**Schwere:** Gering (sprachlich)
**Ort:** EN, Remark "Odd $n$", Abschnitt 3, Zeile ~234–238
**Befund:** Der Text enthält einen unvollendeten Satz:
> "giving $n = 2 + 2 = 4$, odd only for... actually this shows directly..."
Das "odd only for..." ist abgebrochen und wird durch "actually this shows..." korrigiert. Der erste Halbsatz ist sinnlos stehen gelassen.
**Korrektur:** Streichen: "giving $n = 2 + 2 = 4$, odd only for..." (gesamter Halbsatz), da der folgende Satz korrekt erklärt.

#### ⚠️ BUG-B4-P20-EN-002 — Systematisch falsche Zahlenwerte in Tabelle (EN)
**Schwere:** Hoch (Rechenfehler, mehrere Einträge betroffen)
**Ort:** EN, Abschnitt 4 (Numerical Values), Tabelle Zeilen ~247–258
**Befund:** Die Tabelle enthält für **jeden Eintrag** einen zusätzlichen Faktor $(p-1)/(p-2)$ für eine Primzahl $p$, die $n$ **nicht** teilt. Korrekte Formel: $\mathfrak{S}(n) = 2C_2 \prod_{p \mid n,\, p \geq 3} (p-1)/(p-2)$.

| $n$ | Tabelleninhalt | Korrekt |
|-----|----------------|---------|
| 4 | $2C_2 \cdot (2/1) = 4C_2 \approx 2.641$ | $2C_2 \approx 1.320$ (keine unger. Primteiler) |
| 6 | $2C_2 \cdot (2/1) \cdot (4/3) = 16C_2/3 \approx 3.521$ | $2C_2 \cdot 2/1 = 4C_2 \approx 2.641$ |
| 8 | $2C_2 \cdot (2/1) = 4C_2 \approx 2.641$ | $2C_2 \approx 1.320$ (8 = 2³) |
| 10 | $2C_2 \cdot (2/1) \cdot (4/3) = 16C_2/3 \approx 3.521$ | $2C_2 \cdot 4/3 = 8C_2/3 \approx 1.760$ |
| 12 | $2C_2 \cdot (2/1)\cdot(4/3)\cdot(6/5) \approx 4.225$ | $2C_2 \cdot 2/1 = 4C_2 \approx 2.641$ |
| 30 | $2C_2 \cdot (2/1)\cdot(4/3)\cdot(6/5) \approx 5.1$ | $2C_2 \cdot (2/1)\cdot(4/3) = 16C_2/3 \approx 3.521$ |

Das Muster: Die Tabelle enthält systematisch eine Primzahl mehr als $n$ hat, d.h. es wird immer der nächste nicht teilende Primfaktor $(p_\text{next}-1)/(p_\text{next}-2)$ fälschlich eingeschlossen.
**Korrektur:** Alle Tabelleneinträge neu berechnen anhand $2C_2 \prod_{p \mid n,\, p \geq 3} (p-1)/(p-2)$.

#### ⚠️ BUG-B4-P20-DE-001 — Identische Tabellenfehler (DE)
**Schwere:** Hoch (gleiche Fehler wie EN)
**Ort:** DE, Abschnitt 4, Tabelle
**Befund:** Identische falsche Werte wie in der EN-Tabelle.
**Korrektur:** Analog zu EN.

#### Numerische Tabelle im Goldbach-Kometen-Abschnitt (Abschnitt 7 EN / Abschnitt 8 DE)
Die dort angegebenen Prognose-Werte für $r_2(n)$ sind plausibel (Größenordnung korrekt, da es sich um gerundete Schätzwerte handelt). Keine weiteren Fehler identifiziert.

#### GRH-Satz (Theorem 6.2 EN / Satz 7.1 DE)
**EN:** Der Beweis enthält einen zunächst verwirrenden Schritt: Ein "naiver" GRH-Bound $O(N^{3/2})$ wird als zu groß beschrieben — um dann den korrekten Weg zu erklären (GRH verbessert die Hauptbogen-Approximation, nicht nur die Nebenbogen). Das Ergebnis $R_2(n) = \mathfrak{S}(n) \cdot n + O(\sqrt{n} (\log n)^2)$ ist korrekt. **Nicht-blockierender Stilmangel.**
**DE:** Satz ohne ausführlichen Beweis — stattdessen direktes Statement. Korrekt und vollständig für DE.

### Hardy–Littlewood Hypothesis H
- Admissibilitätsdefinition: korrekt.
- Singulärreihendefinition $\mathfrak{S}(\mathcal{H}) = \prod_p (1-1/p)^{-k}(1-\nu_p/p)$: korrekt.
- Spezialfälle (Primzahlsatz, Zwillingsprimes, Cousin-Primes, Tripel): alle korrekt.
- Verknüpfung mit Goldbach (Tupel $\{0, n\}$): korrekt.

### Vollständigkeit
**EN:** Sehr vollständig. Alle Abschnitte vorhanden: Herleitung $\mathfrak{S}(n)$, Euler-Produkt, Positivität, Numerik, Zwillingsprimzahl-Konstante, Hypothesis H, Paritätsproblem, GRH-Ergebnis, Goldbach-Komet, Chen-Bezug.
**DE:** Äquivalent zu EN. Alle Abschnitte vorhanden. Kein wesentlicher Inhalt fehlt.

### Urteil Paper 20
**EN:** ⚠️ ÜBERARBEITUNG NOTWENDIG — BUG-B4-P20-EN-002 (Tabellenfehler, mehrere Einträge), BUG-B4-P20-EN-001 (Satzfragment).
**DE:** ⚠️ ÜBERARBEITUNG NOTWENDIG — BUG-B4-P20-DE-001 (Tabellenfehler).

---

## Zusammenfassung: Gesamtbewertung aller 4 Batches

### Batch 1 (Papers 1–7): Giuga/Lehmer-Serie
**Letzter Review:** Build 5 (2026-03-11)

| Paper | EN | DE | Offene Fehler |
|-------|----|----|---------------|
| 1 — Giuga 3-Prim | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 2 — Lehmer 3-Prim | ✅ DRUCKREIF | ⚠️ BUG-010 | Strukturelle Redundanz (stilistisch) |
| 3 — Giuga-Carmichael | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 4 — Giuga quadratfrei | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 5 — Giuga kein Semiprim | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 6 — Lehmer quadratfrei | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 7 — Lehmer kein Semiprim | ✅ DRUCKREIF | ⚠️ BUG-011 | Schwacher Zwischenschritt (stilistisch) |

### Batch 2 (Papers 8–12): Wilson-Theorem-Serie
**Letzter Review:** Build 5 (2026-03-11)

| Paper | EN | DE | Offene Fehler |
|-------|----|----|---------------|
| 8 — Wilson Grundsatz | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 9 — Wilson Primzahlpotenzen | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 10 — Wilson Quotient | ✅ DRUCKREIF | ✅ DRUCKREIF | BUG-P10-DE-LATEX behoben |
| 11 — Wilson abelsche Gruppen | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 12 — Wilson Anwendungen | ✅ DRUCKREIF | ✅ DRUCKREIF | – |

### Batch 3 (Papers 13–16): Sieblehre-Serie
**Letzter Review:** Build 5 (2026-03-11)

| Paper | EN | DE | Offene Fehler |
|-------|----|----|---------------|
| 13 — Bruns Satz | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 14 — Selbergs Sieb | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 15 — Große Sieb-Ungleichung | ✅ DRUCKREIF | ✅ DRUCKREIF | Build-Fix behoben |
| 16 — Chens Satz | ✅ DRUCKREIF | ✅ DRUCKREIF | Bibitem+Build-Fix behoben |

### Batch 4 (Papers 17–20): Kreismethode-Serie
**Dieser Review:** Build 6 (2026-03-11)

| Paper | EN | DE | Offene Fehler |
|-------|----|----|---------------|
| 17 — Kreismethode | ✅ DRUCKREIF | ✅ DRUCKREIF | – |
| 18 — Vinogradov 3-Prim | ⚠️ | ⚠️ | BUG-B4-P18-EN-001, DE-001, DE-002 |
| 19 — Waringsches Problem | ⚠️ | ⚠️ | BUG-B4-P19-EN-001, DE-001, -002, -003 |
| 20 — Goldbach sing. Reihe | ⚠️ | ⚠️ | BUG-B4-P20-EN-001, -002, DE-001 |

---

## Vollständige Bug-Liste Batch 4 (neue Fehler)

| ID | Paper | Schwere | Beschreibung |
|----|-------|---------|--------------|
| BUG-B4-P18-EN-001 | 18 EN | Mittel | Falsches Vorzeichen in Remark 4.3: `$1+c_2(n)/(2-1)^3$` statt `$1-c_2(n)/(2-1)^3$` |
| BUG-B4-P18-DE-001 | 18 DE | Mittel | Gleiches Vorzeichenproblem in Bemerkung |
| BUG-B4-P18-DE-002 | 18 DE | Gering | Beweis Satz 5.1 (sing. Integral): nur $n \leq N$ Fall, fehlen $N \leq n \leq 2N$ usw. |
| BUG-B4-P19-EN-001 | 19 EN | Gering | Falsche Querverweise: "Paper 17, Def. 5.1" (→ 5.2) und "Paper 17, Lemma 5.3" (→ 5.5) |
| BUG-B4-P19-DE-001 | 19 DE | Mittel | Fehlendes Korollar "Weyl-Schranke G(k)" mit Beweis |
| BUG-B4-P19-DE-002 | 19 DE | Gering | Korollar (Wooley-Schranke) ohne Beweis |
| BUG-B4-P19-DE-003 | 19 DE | Mittel | Kein formaler Satz für effizientes Kongruenzrechnen (nur Bemerkungen) |
| BUG-B4-P20-EN-001 | 20 EN | Gering | Unvollständiger Satz in Remark "Odd n" |
| BUG-B4-P20-EN-002 | 20 EN | Hoch | Numerische Tabelle $\mathfrak{S}(n)$: alle Werte systematisch falsch (1 Extra-Faktor pro Zeile) |
| BUG-B4-P20-DE-001 | 20 DE | Hoch | Gleiche Tabellenfehler wie EN |

---

## Handlungsempfehlungen (Priorität)

### Priorität 1 (Mathematische Fehler)
1. **Paper 20 EN+DE — Tabelle korrigieren** (BUG-B4-P20-EN-002, DE-001): Korrekte Werte:
   - $n=4$: $2C_2$; $n=6$: $4C_2$; $n=8$: $2C_2$; $n=10$: $8C_2/3$; $n=12$: $4C_2$; $n=30$: $16C_2/3$
2. **Paper 18 EN+DE — Vorzeichen korrigieren** (BUG-B4-P18-EN/DE-001): `$+c_2(n)$` → `$-c_2(n)$`

### Priorität 2 (Vollständigkeit)
3. **Paper 18 DE — Beweis Satz 5.1 erweitern** (BUG-B4-P18-DE-002)
4. **Paper 19 DE — Theorem 6.1 (Efficient Congruencing) ergänzen** (BUG-B4-P19-DE-003)
5. **Paper 19 DE — Korollar 4.1 mit Beweis ergänzen** (BUG-B4-P19-DE-001)

### Priorität 3 (Kleinigkeiten)
6. **Paper 20 EN — Satzfragment entfernen** (BUG-B4-P20-EN-001)
7. **Paper 19 EN — Querverweise korrigieren** (BUG-B4-P19-EN-001)
8. **Paper 19 DE — Korollar-Beweis ergänzen** (BUG-B4-P19-DE-002)

---

*Gutachter: Claude (claude-sonnet-4-6) | Datum: 2026-03-11 | Build 6*
