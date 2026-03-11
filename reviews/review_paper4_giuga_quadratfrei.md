# Gutachten: Paper 4 — Jedes Giuga-Pseudoprim ist quadratfrei

**Papers:** `paper4_giuga_squarefree.tex` (EN) · `paper4_giuga_quadratfrei_de.tex` (DE)
**Gutachter:** Claude Code (automatisiertes Peer-Review)
**Datum:** 2026-03-11
**Gesamturteil:** ✅ Annahme empfohlen (beide Versionen)

---

## 1. Zusammenfassung

Beide Papers beweisen, dass jede zusammengesetzte Zahl, die nur die **schwache Giuga-Bedingung** (W)/(S) für alle Primteiler erfüllt, notwendigerweise **quadratfrei** ist. Der Beweis ist elementar und vollständig: Aus $p^2 \mid n$ folgt $p \mid \frac{n}{p}$, damit $\frac{n}{p} - 1 \equiv -1 \pmod{p}$, und (W) erzwingt $p \mid -1$ — Widerspruch, da $p \ge 2$.

---

## 2. Mathematische Korrektheit

### Hauptsatz (Theorem 2.1 / Satz 2.1)
**Korrekt.** Der Beweis per Widerspruch ist lückenlos:

1. Annahme: $p^2 \mid n$ für eine Primzahl $p$. ✓
2. Dann $p \mid \frac{n}{p}$ (da $p^2 \mid n \Rightarrow p \mid \frac{n}{p}$). ✓
3. Also $\frac{n}{p} \equiv 0 \pmod{p}$, d.h. $\frac{n}{p} - 1 \equiv -1 \pmod{p}$. ✓
4. Schwache Bedingung (W): $p \mid \frac{n}{p} - 1$, also $p \mid -1$. ✓
5. Widerspruch: $p \ge 2 \Rightarrow p \nmid -1$. ✓

### Korollare
Beide Korollare folgen direkt aus dem Satz. ✓

### Numerische Illustration
Alle Faktorisierungen wurden nachgeprüft:
- $30 = 2 \cdot 3 \cdot 5$ ✓
- $858 = 2 \cdot 3 \cdot 11 \cdot 13$ ✓ ($6 \cdot 143 = 858$)
- $1722 = 2 \cdot 3 \cdot 7 \cdot 41$ ✓ ($6 \cdot 287 = 1722$)
- $66198 = 2 \cdot 3 \cdot 11 \cdot 17 \cdot 59$ ✓ ($1122 \cdot 59 = 66198$)

Die Aussage, dass keine dieser Zahlen ein Giuga-*Pseudoprim* ist (jede verletzt die starke Bedingung für mind. einen Primfaktor), ist korrekt — es handelt sich um *Giuga-Zahlen* (nur (W)).

---

## 3. Zitate und Literatur

| Referenz | Bewertung |
|---|---|
| Giuga 1950 | ✅ Korrekte Originalreferenz |
| Borwein, Borwein, Girgensohn & Parnes 1996 | ✅ Korrekt — **BUG-001 aus Paper 1 ist hier behoben!** |
| Ingwer 2026a, 2026c | Interne Cross-Referenzen, korrekt |

---

## 4. Vergleich EN ↔ DE

| Aspekt | EN | DE |
|---|---|---|
| Mathematik | Identisch ✓ | Identisch ✓ |
| Bezeichnungen | (W) = schwach, (S) = stark | (S) = schwach, (St) = stark |
| Struktur | Identisch | Identisch |

Die unterschiedlichen Bezeichnungen (W)/(S) im EN vs. (S)/(St) im DE sind konsistent mit den bisherigen Papers und erzeugen keine Verwirrung im jeweiligen Kontext.

---

## 5. Qualität und Darstellung

- **Beweiseleganz:** Sehr hoch. Einzeiliger Kern-Widerspruch.
- **Didaktik:** Die zweite Bemerkung erklärt die Intuition nochmals in Prosa — sinnvoll.
- **Einleitung:** Sauber, gibt Kontext zur Giuga-Vermutung und Borwein et al. (1996).
- **Abschnitt "Consequences":** Die Korollare sind trivial aber nützlich als Referenzpunkte.

---

## 6. Befunde

Keine Fehler gefunden. Beide Versionen sind publishable.

---

## 7. Empfehlung

**Annahme ohne Auflagen** (beide Sprachversionen).
Das Paper liefert das essentielle strukturelle Fundament für alle weiteren Papers in dieser Reihe.
