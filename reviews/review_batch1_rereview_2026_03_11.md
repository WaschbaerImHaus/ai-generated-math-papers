# Re-Review-Bericht: Batch 1 — Giuga/Lehmer-Serie (Papers 1–7)

**Datum:** 2026-03-11
**Gutachter:** Claude Sonnet 4.6 (Re-Review)
**Build:** 5
**Basis:** Dateien in `papers/reviewed/batch1/`

---

## Überblick

Alle 14 Dateien (Papers 1–7, je EN+DE) wurden vollständig neu gelesen und auf mathematische Korrektheit, Vollständigkeit, Stimmigkeit und Beweis-Eigenschaften überprüft.

**Gesamturteil:** Die Hauptaussagen aller Papers sind mathematisch korrekt. Zwei Nebenbefunde wurden identifiziert.

---

## Paper 1 — Kein Giuga-Pseudoprim mit 3 Primfaktoren
**Dateien:** `paper1_giuga_three_primes.tex` (EN), `paper1_giuga_drei_primfaktoren_de.tex` (DE)

### Prüfung

**Lemma Quadratfreiheit:**
- Wenn $p^2 \mid n$, dann $p \mid (n/p)$, also $(n/p) - 1 \equiv -1 \pmod{p}$, daher $p \mid (-1)$ — Widerspruch. ✅

**Lemma kein 2-Primfaktor-Produkt:**
- Aus $q \mid (p-1)$ und $p < q$ folgt $p - 1 = 0$, also $p = 1$ — kein Prim. ✅

**Fall $n = 2qr$ (Paritäts-Hindernis):**
- Aus $r = 2q - 1$ folgt $r - 1 = 2(q-1)$ — gerade. Aber $2q - 1$ ist ungerade → Widerspruch. ✅

**Fall $n = pqr$ (alle ungerade):**
- $\text{lcm}(r, r-1) = r(r-1) \geq (q+2)(q+1) > pq - 1$. Vollständige Abschätzung korrekt. ✅

**Konsistenz EN ↔ DE:** Notation (W)/(S) vs. (S)/(St) unterschiedlich — stilistisch, mathematisch äquivalent. ✅

**LaTeX:** Alle Umgebungen definiert, alle Bibitems vorhanden. ✅

### Urteil
**DRUCKREIF** — keine offenen Fehler. (BUG-001 aus Vorrevision behoben.)

---

## Paper 2 — Lehmer-Totient für 3-Primfaktor-Produkte
**Dateien:** `paper2_lehmer_three_primes.tex` (EN), `paper2_lehmer_drei_primfaktoren_de.tex` (DE)

### Prüfung

**EN-Version:**
- Proposition Semiprim: $(p-1)(q-1) \mid (p+q-2)$ → für $p \geq 3$: $(p-1)(q-1) > p+q-2$. Beweis vollständig. ✅
- Theorem $n = 2qr$: Fallunterscheidung $b = 1$ ($r = 2q$ zusammengesetzt) und $b \geq 2$ ($q \geq r$ Widerspruch). ✅
- Theorem $n = pqr$: Schlüsselidentität und Schrankenargument vollständig. ✅

**DE-Version — BUG-010 (Neubefund):**
- Zeilen 140–142: Vollständiger Paritätsbeweis für den Fall $n = 2qr$ ($(q-1)(r-1)$ gerade, $2qr-1$ ungerade → keine Teilbarkeit).
- Zeilen 144–155: Ein zweiter, alternativer Fallunterscheidungs-Beweis (Fall $b=1$, $b\geq 2$) wird ohne Überleitung angehängt.
- **Befund:** Der abgebrochene `\ldots`-Entwurfsüberrest (laut Vorrevision) ist **nicht mehr vorhanden** — beide Beweise sind nun vollständig. Das Problem ist jetzt nur noch **strukturell**: zwei vollständige Beweise in einem `\begin{proof}`-Block ohne klare "Alternativ:"-Markierung. Dies ist stilistisch unordentlich, mathematisch aber korrekt.
- **Empfehlung:** Entweder Paritätsargument (Zeilen 140–142) als einzigen Beweis behalten, oder alternativen Block mit `\medskip\noindent\textit{Alternativer Beweis:}` einleiten.
- **Status BUG-010:** Der ursprüngliche `\ldots`-Fehler ist behoben. Neue Einstufung: stilistische Redundanz, keine mathematische Blockierung.

**Konsistenz EN ↔ DE:** Mathematisch äquivalent. ✅

**LaTeX:** ✅

### Urteil
**DRUCKREIF** (mit Stilanmerkung) — BUG-010 inhaltlich behoben, strukturelle Bereinigung empfohlen.

---

## Paper 3 — Giuga-Pseudoprime und Carmichael-Zahlen unverträglich
**Dateien:** `paper3_giuga_carmichael.tex` (EN), `paper3_giuga_carmichael_de.tex` (DE)

### Prüfung

**Prop. Carmichael-Zahlen sind ungerade:**
- $(-1)^{n-1} = 1$ für gerades $n$ mit ungeradem $n-1$ erzwingt $-1 \equiv 1 \pmod{n}$, also $n \mid 2$. ✅

**Theorem notwendige Kongruenz:**
- Aus Korselt-Bedingung und Giuga-Bedingung via CRT: $n \equiv p \pmod{p^2(p-1)}$. Beweis korrekt und vollständig. ✅

**Theorem CRT-Hindernis {3,5,7}:**
- Schritt 1: CRT liefert $n \equiv 705 \pmod{900}$. Nachgerechnet: korrekt.
- Schritt 2: $\gcd(900, 294) = 6$; $7 - 705 = -698 \equiv 4 \pmod{6}$, also kein CRT-Lift möglich. Verifiziert: $698 = 116 \cdot 6 + 2$, also $-698 \equiv -2 \equiv 4 \pmod{6}$. ✅

**Vollständigkeit:** Nur {3,5,7} behandelt; {3,5,11} usw. als offene Probleme korrekt ausgewiesen. ✅

**LaTeX:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 4 — Giuga-Pseudoprime sind quadratfrei
**Dateien:** `paper4_giuga_squarefree.tex` (EN), `paper4_giuga_quadratfrei_de.tex` (DE)

### Prüfung

**Satz Quadratfreiheit:** $p^2 \mid n$ → $p \mid (n/p)$ → $(n/p) \equiv 1 \pmod{p}$ → $p \mid (-1)$ — Widerspruch. ✅

**Numerische Illustration:** Bekannte Giuga-Zahlen {30, 858, 1722, 66198} korrekt. ✅

**LaTeX:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 5 — Kein Giuga-Semiprim
**Dateien:** `paper5_giuga_no_semiprime.tex` (EN), `paper5_giuga_kein_semiprim_de.tex` (DE)

### Prüfung

**Theorem:** Kein $n = pq$ (Semiprim) ist Giuga-Pseudoprim.
- Aus $q \mid (p-1)$ und $p - 1 < q$ folgt $p - 1 = 0$ — kein Prim. ✅

**LaTeX:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen. (BUG-005, 006, 007 aus Vorrevision behoben.)

---

## Paper 6 — Lehmer-Zahlen sind quadratfrei
**Dateien:** `paper6_lehmer_squarefree.tex` (EN), `paper6_lehmer_quadratfrei_de.tex` (DE)

### Prüfung

**Satz:** $p^2 \mid n$ → $p \mid \varphi(n)$ → aus $\varphi(n) \mid (n-1)$ folgt $p \mid (n-1)$ → aber $p \mid n$ → $p \mid 1$ — Widerspruch. ✅

**Vergleich mit Giuga-Fall (Paper 4):** Korrekter Hinweis auf strukturelle Ähnlichkeit. ✅

**LaTeX:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 7 — Kein Lehmer-Semiprim
**Dateien:** `paper7_lehmer_no_semiprime.tex` (EN), `paper7_lehmer_kein_semiprim_de.tex` (DE)

### Prüfung

**Schlüsselidentität:** $pq - 1 = (p-1)(q-1) + (p-1) + (q-1)$. Korrekt. ✅

**Theorem (Fall $p = 2$):** $(q-1) \mid q$ → $(q-1) \mid 1$ → $q = 2 = p$ — Widerspruch. ✅

**Theorem (Fall $p \geq 3$):**
- EN-Version (Zeile 102): $(p-2)(q-2) - 1 \geq 1 \cdot 2 - 1 = 1 > 0$ — klar und korrekt (nutzt $p-2 \geq 1$ und $q-2 \geq 3$). ✅
- DE-Version (Zeile 98): $(p-2)(q-2) - 1 \geq 1 \cdot (p-1) - 1 = p-2 \geq 1$

**Neuer Befund BUG-011 (DE, Zeile 98):**
Der Zwischenschritt $1 \cdot (p-1)$ impliziert $q - 2 \geq p - 1$. Da $q \geq p+2$ (zuvor bewiesen), gilt $q-2 \geq p > p-1$ — die Aussage ist korrekt, der Schritt aber unnötig schwach und inkonsistent mit der EN-Version. **Empfehlung:** Ersetzen durch $\geq 1 \cdot p - 1 = p - 1 \geq 2 > 0$ (analog zur EN-Version mit konkreten Schranken).

**Konsistenz EN ↔ DE:** Mathematisch korrekt, aber Zwischenschritt in DE weniger präzise. ✅ (mit Anmerkung)

**LaTeX:** ✅

### Urteil
**DRUCKREIF** (mit Präzisierungsempfehlung) — BUG-011 dokumentiert, kein mathematischer Fehler.

---

## Zusammenfassung Batch 1 Re-Review

| Paper | Titel | EN | DE | Urteil |
|---|---|---|---|---|
| 1 | Giuga 3-Prim | ✅ | ✅ | DRUCKREIF |
| 2 | Lehmer 3-Prim | ✅ | ✅ | DRUCKREIF (Stilanmerkung) |
| 3 | Giuga-Carmichael | ✅ | ✅ | DRUCKREIF |
| 4 | Giuga quadratfrei | ✅ | ✅ | DRUCKREIF |
| 5 | Giuga kein Semiprim | ✅ | ✅ | DRUCKREIF |
| 6 | Lehmer quadratfrei | ✅ | ✅ | DRUCKREIF |
| 7 | Lehmer kein Semiprim | ✅ | ✅ | DRUCKREIF (Stilanmerkung) |

### Befunde dieser Re-Review-Runde
- **BUG-010:** Ursprünglicher `\ldots`-Entwurfsüberrest behoben; strukturelle Redundanz verbleibt (zwei vollständige Beweise ohne klare Trennung).
- **BUG-011 (neu):** Paper 7 DE, Zeile 98 — Zwischenschritt der Ungleichungskette schwächer als nötig, inhaltlich aber korrekt.

### Alle Hauptsätze mathematisch korrekt ✅
