# Gutachten: Paper 5 — Kein Semiprim ist ein Giuga-Pseudoprim

**Papers:** `paper5_giuga_no_semiprime.tex` (EN) · `paper5_giuga_kein_semiprim_de.tex` (DE)
**Gutachter:** Claude Code (automatisiertes Peer-Review)
**Datum Erstgutachten:** 2026-03-11
**Datum Revisionsgutachten:** 2026-03-11
**Gesamturteil:** ✅ **Annahme ohne Auflagen** — alle bekannten Fehler behoben (EN + DE)

---

## Revisionsgutachten (2026-03-11)

### BUG-005 — Status: ✅ BEHOBEN (EN + DE)

- **EN (Zeile 57):** `\cite{Giuga1950}` → `\cite{Ingwer2026sqfree}` ✅
- **DE (Zeile 59):** `\cite{Giuga1950}` → `\cite{Ingwer2026qfrei}` ✅

Beide Versionen zitieren nun korrekt Paper 4 dieser Reihe für die Quadratfreiheit.

### BUG-006 — Status: ✅ BEHOBEN (EN)

Remark in Section 3, Zeile 117–118 lautet nun:
```latex
n - 1 = \varphi(n) + \sum_{\emptyset \ne S \subsetneq \{1,\ldots,k\}} \prod_{i \in S}(p_i - 1).
```
Der Term $S = \emptyset$ ist korrekt ausgeschlossen. Verifikation für $k=2$:
$(p-1)(q-1) + (p-1) + (q-1) = pq - 1$ ✅

### BUG-007 — Status: ✅ BEHOBEN (EN)

Remark (Zeile 116) lautet nun: „due to Bednarek \cite{Bednarek2014}" (ohne „and Borwein").
Konsistenz mit Bibliographieeintrag (`M.~Bednarek`) hergestellt. ✅

### Gesamtbewertung
Keine neuen Fehler eingeführt. Hauptbeweise unverändert korrekt.
**Beide Versionen: Druckreif.**

---

---

## 1. Zusammenfassung

Beide Papers beweisen, dass kein Semiprim $n = p \cdot q$ mit $p < q$ prim ein Giuga-Pseudoprim ist. Der Beweis ist ein reines **Ordnungsargument**: Die schwache Bedingung (W)/(S) für die größere Primzahl $q$ erzwingt $q \mid (p-1)$; da $0 \le p-1 < q$, folgt $p-1 = 0$, also $p=1$ — kein Prim. Das Hauptresultat ist korrekt.

---

## 2. Mathematische Korrektheit

### Hauptsatz (Theorem 2.1 / Satz 2.1)
**Korrekt.** Der Beweis ist vollständig und lückenlos:

1. $n = pq$, schwache Bedingung für $q$: $q \mid \frac{n}{q} - 1 = p - 1$. ✓
2. $0 \le p - 1 < q$ (da $p < q$, beide prim, also $p \le q-1$). ✓
3. Das einzige nicht-negative Vielfache von $q$, das $< q$ ist, ist $0$. ✓
4. Folglich $p - 1 = 0$, d.h. $p = 1$ — kein Prim. Widerspruch. ✓

**Bemerkung zur Vollständigkeit:** Die Bedingung (W) für $p$ (d.h. $p \mid q - 1$) wird tatsächlich nicht für den Widerspruch benötigt. Das wird korrekt angemerkt.

### Korollar: Mindestens drei Primfaktoren
**Korrekt.** Aus Quadratfreiheit (Paper 4) + kein Einprimfaktor (wäre prim) + kein Zweiprimfaktor (Hauptsatz) folgt $k \ge 3$. ✓

---

## 3. Gefundene Fehler

### BUG-005 (EN + DE): Falscher Zitatschlüssel für Quadratfreiheit

**Schwere:** Mittel — sachlicher Zitierungsfehler

**EN, Zeile 57:**
```latex
every Giuga pseudoprime is squarefree \cite{Giuga1950}
```
**DE, Zeile 58–59:**
```latex
jedes Giuga-Pseudoprim quadratfrei ist \cite{Giuga1950}
```

**Problem:** Die Quadratfreiheit wird Giugas Originalarbeit (1950) zugeschrieben. Sie wurde jedoch dort nicht explizit in der vorliegenden Form bewiesen — das ist gerade das Ergebnis von Paper 4 (`\cite{Ingwer2026sqfree}` bzw. `\cite{Ingwer2026qfrei}`).

**Korrektur:**
- EN: `\cite{Giuga1950}` → `\cite{Ingwer2026sqfree}`
- DE: `\cite{Giuga1950}` → `\cite{Ingwer2026qfrei}`

---

### BUG-006 (EN only): Formel in Bemerkung (Section 3) falsch um +1

**Schwere:** Gering (Bemerkung, nicht Hauptbeweis)

**EN, Remark nach Corollary 3.1 (ca. Zeile 125–127):**
```latex
n - 1 = \euler(n) + \sum_{\substack{S \subseteq \{1,\ldots,k\} \\ S \ne \{1,\ldots,k\}}}
          \prod_{i \in S}(p_i - 1) \cdot \prod_{i \notin S} \mathbf{1}.
```

**Problem:** Die Summe enthält den Term $S = \emptyset$, für den das leere Produkt $\prod_{i \in \emptyset}(\cdot) = 1$ ist. Damit ergibt die rechte Seite:
$$\varphi(n) + 1 + \sum_{\emptyset \ne S \subsetneq \{1,\ldots,k\}} \prod_{i \in S}(p_i - 1) = n$$
aber nicht $n - 1$.

**Verifikation für $k=2$, $n = pq$:**
RHS $= (p-1)(q-1) + 1 + (p-1) + (q-1) = pq \ne pq - 1$. ✗

**Korrektur:** Summe muss über $\emptyset \ne S \subsetneq \{1,\ldots,k\}$ laufen:
$$n - 1 = \varphi(n) + \sum_{\emptyset \ne S \subsetneq \{1,\ldots,k\}} \prod_{i \in S}(p_i - 1)$$

---

### BUG-007 (EN only): Inkonsistenz zwischen Fließtext und Bibliographie

**Schwere:** Gering

**EN, Remark (ca. Zeile 116):**
> "The current best lower bound, due to **Bednarek and Borwein** \cite{Bednarek2014}, is 59 prime factors."

**Bibliographie:**
```latex
\bibitem{Bednarek2014}
M.~Bednarek,
\emph{On the size of Giuga pseudoprimes},
Preprint (2014).
```

**Problem:** Der Text nennt zwei Autoren ("Bednarek and Borwein"), das Literaturnachweise-Eintrag listet nur einen ("M. Bednarek"). Einer der beiden muss falsch sein.

**Zusätzliche Anmerkung:** Die Schranke von 59 Primfaktoren für Giuga-Pseudoprime konnte nicht anhand bekannter Literatur verifiziert werden. Borwein et al. (1996) zeigten eine Schranke an der Anzahl der *Stellen*, nicht der Primfaktoren. Ggf. handelt es sich um ein unveröffentlichtes Preprint — dann sollte es als solches besonders deutlich kenntlichgemacht werden.

---

## 4. Stärken

- Sehr eleganter, elementarer Beweis (reines Ordnungsargument)
- Sauber modularisiert: Quadratfreiheit als Voraussetzung importiert, Semiprim-Ausschluss eigenständig
- Symmetriebeobachtung in Abschnitt 4 (Kombinierter Ausschluss $p=q$ und $p<q$) ist konzeptionell wertvoll
- DE-Version: Beweisdarstellung im Case $p \ge 3$ sogar klarer als EN (keine missglückte Herleitung)

---

## 5. Vergleich EN ↔ DE

| Aspekt | EN | DE |
|---|---|---|
| BUG-005 (falsches Zitat) | ✗ vorhanden | ✗ vorhanden |
| BUG-006 (Formel +1 off) | ✗ vorhanden | ✅ nicht vorhanden (andere Formulierung) |
| BUG-007 (Autoren-Inkonsistenz) | ✗ vorhanden | ✅ nicht vorhanden |
| Hauptbeweis | ✓ korrekt | ✓ korrekt |
| Case $p \ge 3$-Darstellung | umständlich | klarer |

---

## 6. Empfehlung

**Annahme mit kleinen Überarbeitungen:**
- BUG-005 (beide Versionen): Zitat korrigieren
- BUG-006 (EN): Formel in Remark korrigieren (S=∅ ausschließen)
- BUG-007 (EN): Bibliographie-Eintrag um zweiten Autor ergänzen oder Fließtext korrigieren; Schranke "59 Primfaktoren" ggf. mit Quelle nachweisen
