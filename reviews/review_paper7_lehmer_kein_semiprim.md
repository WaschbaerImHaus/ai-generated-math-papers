# Gutachten: Paper 7 — Kein Semiprim erfüllt Lehmers Totient-Bedingung

**Papers:** `paper7_lehmer_no_semiprime.tex` (EN) · `paper7_lehmer_kein_semiprim_de.tex` (DE)
**Gutachter:** Claude Code (automatisiertes Peer-Review)
**Datum Erstgutachten:** 2026-03-11
**Datum Revisionsgutachten:** 2026-03-11
**Gesamturteil:** ✅ **Annahme ohne Auflagen** — alle bekannten Fehler behoben (EN + DE)

---

## Revisionsgutachten (2026-03-11)

### BUG-008 — Status: ✅ BEHOBEN (EN)

Remark in Section 3, Zeile 117–118 lautet nun:
```latex
n - 1 = \varphi(n) + \sum_{\emptyset \ne S \subsetneq \{1,\ldots,k\}} \prod_{i \in S}(p_i - 1).
```
Korrekt: nur nicht-leere echte Teilmengen. Verifikation für $k=2$:
$(p-1)(q-1) + (p-1) + (q-1) = pq - 1$ ✅

### BUG-009 — Status: ✅ BEHOBEN (EN)

Der Redaktionsüberrest „Hmm — let us handle $p=2$ separately." wurde vollständig entfernt.
Der eingebettete Beweis `proof of (eq:ineq)` behandelt die Fälle $p=2$ und $p \ge 3$
nun sauber als nummerierte Fälle ohne informellen Kommentar. ✅

### Gesamtbewertung
Keine neuen Fehler eingeführt. Hauptbeweise unverändert korrekt.
**Beide Versionen: Druckreif.**

---

---

## 1. Zusammenfassung

Beide Papers beweisen, dass kein Semiprim $n = p \cdot q$ mit $p < q$ prim die Lehmer-Bedingung $\varphi(n) \mid (n-1)$ erfüllt. Die zentrale Idee ist die Schlüsselidentität:
$$pq - 1 = (p-1)(q-1) + (p-1) + (q-1)$$
Dies zeigt $pq - 1 \equiv p + q - 2 \pmod{(p-1)(q-1)}$, und man weist nach, dass $(p-1)(q-1) \nmid (p+q-2)$ durch Fallunterscheidung $p=2$ vs. $p \ge 3$.

---

## 2. Mathematische Korrektheit

### Schlüsselidentität (Gleichung eq:ident)
**Korrekt.**
$$(p-1)(q-1) + (p-1) + (q-1) = pq - p - q + 1 + p - 1 + q - 1 = pq - 1 \checkmark$$

### Hauptsatz — Fall $p = 2$
**Korrekt.**
Bedingung: $(q-1) \mid (2q-1)$.
Da $2q - 1 = 2(q-1) + 1$, folgt $(q-1) \mid 1$, also $q = 2 = p$. Widerspruch zu $q > p$. ✓

### Hauptsatz — Fall $p \ge 3$

**EN-Version (Zeilen 99–112):**

Berechnung: $(p-1)(q-1) - (p+q-2) = pq - 2p - 2q + 3$.

Faktorisierung: $pq - 2p - 2q + 3 = (p-2)(q-2) - 1$. ✓

Für $p \ge 3$ (beide Primzahlen ungerade, $q > p$, also $q \ge p+2$):
- $p - 2 \ge 1$ ✓
- $q - 2 \ge q - p \ge 2$ ✓ (da $q \ge p+2$)
- $(p-2)(q-2) \ge 1 \cdot 2 = 2$, also $(p-2)(q-2) - 1 \ge 1 > 0$ ✓

Folglich $(p-1)(q-1) > p + q - 2$, also $(p-1)(q-1) \nmid (p+q-2)$. ✓

**DE-Version (Zeile 94–99):**

Schärfere Abschätzung: $q \ge p+2 \Rightarrow q-2 \ge p$, also:
$(p-2)(q-2) \ge 1 \cdot (p-1) = p-1$, woraus $(p-2)(q-2) - 1 \ge p - 2 \ge 1$. ✓

*Anmerkung:* Die DE-Version nutzt $q-2 \ge p-1$ (schwächer als das bewiesene $q-2 \ge p$), was die Rechnung leicht unscharf macht, aber korrekt bleibt.

### Korollar: Mindestens drei Primfaktoren
**Korrekt.** Aus Paper 6 (Quadratfreiheit) + kein Einprimfaktor + kein Zweiprimfaktor (Hauptsatz) folgt $k \ge 3$. ✓

---

## 3. Gefundener Fehler

### BUG-008 (EN only): Formel in Remark (Section 3) falsch um +1

**Schwere:** Gering (Bemerkung, nicht Hauptbeweis)

**EN, Remark (ca. Zeilen 125–127):**
```latex
n - 1 = \euler(n) + \sum_{\substack{S \subseteq \{1,\ldots,k\} \\ S \ne \{1,\ldots,k\}}}
          \prod_{i \in S}(p_i - 1) \cdot \prod_{i \notin S} \mathbf{1}.
```

**Problem:** Die Summe enthält den Term $S = \emptyset$, für den das leere Produkt = 1. Damit:
$$\text{RHS} = \varphi(n) + 1 + \sum_{\emptyset \ne S \subsetneq \{1,\ldots,k\}} \prod_{i \in S}(p_i-1) = n$$
statt $n - 1$.

**Verifikation für $k=2$:** $(p-1)(q-1) + 1 + (p-1) + (q-1) = pq \ne pq-1$. ✗

Die korrekte Formel lautet:
$$n - 1 = \varphi(n) + \sum_{\emptyset \ne S \subsetneq \{1,\ldots,k\}} \prod_{i \in S}(p_i - 1)$$

**Korrektur:** Entweder Summenindex auf $\emptyset \ne S \subsetneq \{1,\ldots,k\}$ einschränken, oder letzten Teil des $\prod_{i \notin S} \mathbf{1}$-Terms weglassen und S=∅ explizit ausschließen.

---

## 4. Anmerkungen zur Darstellung

### EN: Geschachtelte `proof`-Umgebungen
Der EN-Text enthält `\begin{proof}[Verification of ...]` und `\begin{proof}[Proof of ineq]` **innerhalb** des Hauptbeweises. Das ist in `amsart` technisch zulässig, aber stilistisch ungewöhnlich. Eine übersichtlichere Variante wäre, die Verifikation der Identität in eine Bemerkung auszulagern.

### EN: Abgebrochene Herleitung in `proof of (eq:ineq)`
Zeilen 97–101 enthalten eine halbfertige Rechnung mit "Hmm — let us handle $p=2$ separately." Dies ist ein editorischer Überrest. Bei einer Zeitschrifteneinreichung sollte dieser Satz entfernt werden.

### DE: Zweite Bemerkung (Äquivalenz zu $(p-2)(q-2) > 1$)
Die Feststellung $(p-1)(q-1) > p+q-2 \Leftrightarrow (p-2)(q-2) > 1$ ist ein nützliches Kriterium. Die Begründung "da $p-2 \ge 1$ und $q-2 \ge p-1 \ge 2$" ist korrekt (da $p \ge 3$, $q \ge p+2$). ✓

---

## 5. Literatur

| Referenz | Bewertung |
|---|---|
| Lehmer 1932 | ✅ Korrekt |
| Cohen & Hagis 1980 | ✅ Korrekt |
| Ingwer 2026 (sqfree/qfrei) | ✅ Korrekt (internes Cross-Referenz) |
| Ingwer 2026 (lehmer3) | Korrekte interne Vorausreferenz |

---

## 6. Vergleich EN ↔ DE

| Aspekt | EN | DE |
|---|---|---|
| BUG-008 (Formel in Remark) | ✗ vorhanden | ✅ nicht vorhanden (andere Formulierung) |
| Redaktioneller Überrest ("Hmm") | ✗ vorhanden | ✅ nicht vorhanden |
| Geschachtelte Beweise | ✗ stilistisch ungünstig | ✅ klarer strukturiert |
| Hauptbeweis Korrektheit | ✓ | ✓ |

Die DE-Version ist in der Darstellung klarer und enthält keine Fehler.

---

## 7. Empfehlung

**Annahme mit kleinen Überarbeitungen (EN):**
1. **BUG-008**: Formel im Remark korrigieren ($S = \emptyset$ ausschließen)
2. Redaktionellen Überrest "Hmm — let us handle $p=2$ separately." entfernen

**Annahme ohne Auflagen (DE).**
