# Review: Paper 7 — No Semiprime Satisfies Lehmer's Totient Condition (EN + DE)

**Dateien:**
- `papers/reviewed/batch1/paper7_lehmer_no_semiprime.tex` (EN)
- `papers/reviewed/batch1/paper7_lehmer_kein_semiprim_de.tex` (DE)

**Reviewer:** Claude (selbst)
**Datum:** 2026-03-11
**Ergebnis:** ✅ Nach Bugfix (BUG-B1-02 + BUG-B1-02-DE) mathematisch korrekt

---

## Überblick

Paper 7 beweist: Kein Semiprim $n = pq$ ($p < q$ prim) erfüllt Lehmers Bedingung
$\varphi(n) \mid n-1$, d.h.\ $(p-1)(q-1) \mid pq-1$.

Schlüsselidentität: $pq - 1 = (p-1)(q-1) + (p-1) + (q-1)$, also
$pq - 1 \equiv p+q-2 \pmod{(p-1)(q-1)}$.

---

## Beweisstruktur

**Fall p=2:** Bedingung wird $(q-1) \mid q$. Da $q = (q-1)+1$, gilt $(q-1) \mid 1$,
also $q=2=p$ — Widerspruch. ✅

**Fall p≥3:** $(p-1)(q-1) - (p+q-2) = (p-2)(q-2) - 1 \ge 1 > 0$,
also $(p-1)(q-1) > p+q-2 > 0$ und kein Teiler. ✅

---

## Behobene Bugs

**BUG-B1-02 (EN):**
- **Ursprünglich falsch**: „Since $2q-1 = 2(q-1)+1$, we get $(q-1) \mid 1$"
- **Problem**: Die Bedingung für $p=2$ lautet $(q-1) \mid q$, **nicht** $(q-1) \mid (2q-1)$.
  Der Ausdruck $2q-1$ tauchte unvermittelt auf.
- **Korrektur**: „Since $q = (q-1)+1$, we have $q \equiv 1 \pmod{q-1}$,
  so $(q-1) \mid 1$" — korrekte Herleitung aus der richtigen Bedingung.

**BUG-B1-02-DE (DE):**
- Selbes Problem in der deutschen Version: Bedingung fälschlich als $(q-1) \mid (2q-1)$
  angegeben (statt $(q-1) \mid q$).
- Korrigiert analog zur englischen Version.

---

## Bugs

| ID | Status | Beschreibung |
|----|--------|-------------|
| BUG-B1-02 | ✅ Behoben | EN: Bedingung $(q-1)\|q$ fälschlich zu $(q-1)\|(2q-1)$ aufgeschrieben; korrekte Herleitung via $q=(q-1)+1$ eingefügt |
| BUG-B1-02-DE | ✅ Behoben | DE: Selber Fehler in der deutschen Version behoben |

---

## Fazit

Paper 7 beweist das Nicht-Semiprim-Lemma für Lehmer korrekt und vollständig.
Der einzige Fehler war eine falsch aufgeschriebene Zwischenbedingung im Fall $p=2$,
die die richtige Schlussfolgerung aus einem falschen Ausdruck zog.
Nach dem Fix ist der Beweis mathematisch rigoros.
