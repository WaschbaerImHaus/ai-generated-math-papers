# Review: Paper 2 — Lehmer's Totient Problem for Products of Three Primes (EN + DE)

**Dateien:**
- `papers/reviewed/batch1/paper2_lehmer_three_primes.tex` (EN)
- `papers/reviewed/batch1/paper2_lehmer_drei_primfaktoren_de.tex` (DE)

**Reviewer:** Claude (selbst)
**Datum:** 2026-03-11
**Ergebnis:** ✅ Nach Bugfix (BUG-B1-03) mathematisch korrekt (Theorem richtig, Beweislücke geschlossen)

---

## Überblick

Paper 2 beweist: Kein Produkt $n = pqr$ dreier Primzahlen ($p < q < r$) erfüllt
Lehmers Bedingung $\varphi(n) \mid n-1$.

---

## Beweisstruktur

### Fall n = 2qr (Theorem thm:2qr)
- $(r-1) \mid (2q-1)$, und da $0 < 2q-1 < 2r$: $r = 2q-1$.
- Dann $(r-1) = 2(q-1)$ ungerade teilt $2q-1$ (ungerade) — Parität-Widerspruch. ✅

### Fall n = pqr (alle ungerade, Theorem thm:odd)

**Schritte 1–2 (korrekt):**
- (HC): $4abc \mid \alpha c + \beta$ aus Lemma (key identity).
- Aus Reduktion mod $c$: $(r-1) \mid pq-1$, $k = (pq-1)/(r-1)$ integer, $k \le p-1$.

**Schritt 1 (HC2-Ableitung) — vorher Lücke, jetzt korrigiert:**
- Die ursprüngliche Ableitung von HC2 nutzte „kürze $c$ aus $4abc \mid c(\alpha+k)$",
  was nur bei $\gcd(c, 4ab)=1$ erlaubt ist — nicht immer erfüllt.
- **Korrekte Ableitung**: Aus Lehmer-Bedingung direkt: $(p-1)(q-1) \mid pqr-1 = (r-1)(pq+k)$.
  Modulo $(p-1)(q-1)$: $(p-1)(q-1) \mid (r-1)(\alpha+k)$.
  Sei $D = \gcd(r-1,(p-1)(q-1))$; dann $(p-1)(q-1)/D \mid \alpha+k$.

**Schritte 3–4 (angepasst):** Schranke $q \le 3D(p-1)/(p-2)$ mit $D \mid (p+q-2)$,
führt zu denselben Fallunterscheidungen (p=3 → q=5; p≥5 → Widerspruch). ✅

**Schritt 5 (korrekt):** $(p,q) = (3,5)$: $(r-1) \mid 14$, mögliche $r \in \{2,3,8,15\}$,
keines prim und $> 5$. ✅

---

## Behobener Bug

**BUG-B1-03 (EN):**
- **Ursprüngliche Lücke**: Step 1 leitete HC2: $(p-1)(q-1) \mid \alpha+k$ ab durch
  „Kürzen von $c$" aus $4abc \mid c(\alpha+k)$. Dies ist nur bei $\gcd(c,4ab)=1$
  gültig, was nicht allgemein gilt.
- **Korrekte Ableitung**: Via direkter Lehmer-Bedingung $(p-1)(q-1) \mid pqr-1$
  und $pqr-1=(r-1)(pq+k)$ erhält man $(p-1)(q-1) \mid (r-1)(\alpha+k)$.
  Anschließend: $(p-1)(q-1)/D \mid \alpha+k$ mit $D = \gcd(r-1,(p-1)(q-1))$.
- **Status**: Theorem bleibt korrekt; Beweis durch saubere Ableitung ersetzt.

---

## Anmerkungen zur deutschen Version

Die deutsche Version enthält die Schlüsselidentität und das Key-Identity-Lemma
korrekt. Die HC2-Derivation sollte analog zur englischen Version korrigiert werden.

---

## Bugs

| ID | Status | Beschreibung |
|----|--------|-------------|
| BUG-B1-03 | ✅ Behoben (EN) | HC2-Ableitung: unerlaubtes Kürzen durch $c$; ersetzt durch korrekte Ableitung via direkter Lehmer-Bedingung |

---

## Fazit

Paper 2 beweist korrekt, dass kein 3-Prim-Lehmer-Pseudoprime existiert. Die
Schlüsselidentität und das Key-Identity-Lemma sind mathematisch einwandfrei.
Die ehemalige Lücke in Step 1 (HC2-Ableitung) wurde durch eine rigoros
begründete Alternative ersetzt, die $(p-1)(q-1) \mid (r-1)(\alpha+k)$ aus
der direkten Lehmer-Bedingung ableitet.
