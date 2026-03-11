# Review: Paper 12 — Anwendungen des Wilson-Satzes (EN + DE)

**Dateien:**
- `papers/batch2/paper12_wilson_applications.tex` (EN)
- `papers/batch2/paper12_wilson_anwendungen_de.tex` (DE)

**Reviewer:** Claude (selbst)
**Datum:** 2026-03-11
**Ergebnis:** ✅ Nach Bugfix (BUG-B2-09) mathematisch korrekt

---

## Überblick

Paper 12 demonstriert vier klassische Anwendungen des Wilson-Satzes:
1. **Primzahltest via Wilson** — theoretisch perfekt, praktisch ineffizient
2. **Quadratische Reste:** $\left(\frac{-1}{p}\right) = (-1)^{(p-1)/2}$ (erstes Ergänzungssatz)
3. **Wolstenholmes Harmonisches Lemma:** $\sum_{j=1}^{p-1} j^{-1} \equiv 0 \pmod p$
4. **Fermats kleiner Satz via Wilson** (Alternativbeweis)

---

## Analyse der Beweise

### Abschnitt 2 — Primzahltest

**Satz 2.1:** $n$ prim $\iff$ $(n-1)! \equiv -1 \pmod n$.

**Beweis für "$(n-1)!\not\equiv -1$" bei zusammengesetztem $n \geq 6$:**
- Fall $n = ab$, $1 < a < b < n$: Beide $a, b$ kommen in $(n-1)!$ vor, also $n \mid (n-1)!$. ✅
- Fall $n = a^2$, $a \geq 3$: $a$ und $2a$ sind beide in $\{1,\ldots,n-1\}$ (da $2a < a^2$ für $a\geq3$) und $a \neq 2a$. Ihr Produkt $2a^2 = 2n$ erscheint in $(n-1)!$, also $n \mid (n-1)!$. ✅ *(Diese Korrektur wurde in Paper 8 ergänzt und ist hier ebenfalls enthalten.)*
- Sonderfälle $n=4$: $(4-1)! = 6 \equiv 2 \pmod 4 \neq -1$. ✅

**Bewertung:** Vollständig und korrekt.

### Abschnitt 3 — Quadratische Reste

**Satz 3.1:** Für Primzahl $p \geq 3$: $\left(\frac{-1}{p}\right) = (-1)^{(p-1)/2}$.

**Beweis:**
1. $(p-1)! = \prod_{j=1}^{m} j \cdot \prod_{j=1}^{m} (p-j) = (-1)^m (m!)^2$ mit $m = (p-1)/2$. ✅
2. Wilson: $(p-1)! \equiv -1 \pmod p$, also $(-1)^m(m!)^2 \equiv -1$. ✅
3. $-1 \equiv (m!)^2 \cdot (-1)^m \pmod p$. Wenn $m$ gerade ($p \equiv 1 \pmod 4$): $-1 \equiv (m!)^2$ → $-1$ ist QR. Wenn $m$ ungerade ($p \equiv 3 \pmod 4$): $(m!)^2 \equiv 1$ → $-1$ ist NQR. ✅

**Bewertung:** Eleganter, vollständiger Beweis. Korrekt. ✅

### Abschnitt 4 — Wolstenholmes Harmonisches Lemma

**Satz 4.1:** Für Primzahl $p \geq 3$: $\displaystyle\sum_{j=1}^{p-1} \frac{1}{j} \equiv 0 \pmod p$.

**Beweis:**
$H = \sum_{j=1}^{p-1} j^{-1}$. Da $j \mapsto j^{-1}$ eine Bijektion auf $\{1,\ldots,p-1\}$ ist, gilt $H \equiv \sum_{j=1}^{p-1} j = p(p-1)/2 \equiv 0 \pmod p$. ✅

*(Alternativ: Paarung $j^{-1} + (p-j)^{-1} \equiv 0 \pmod p$, also $H \equiv 0$.)*

**Bewertung:** Korrekt und elegant. ✅

**Wolstenholmes stärkerer Satz (1862):** Für $p \geq 5$: $H \equiv 0 \pmod{p^2}$. Als Verweis korrekt angegeben. ✅

### Abschnitt 5 — Fermats kleiner Satz via Wilson

**Satz 5.1:** Für Primzahl $p$ und $\gcd(a,p)=1$: $a^{p-1} \equiv 1 \pmod p$.

**Beweis:**
1. $\{a, 2a, \ldots, (p-1)a\}$ ist eine Permutation von $\{1, 2, \ldots, p-1\}$ mod $p$. ✅
2. Produkte: $a^{p-1} \cdot (p-1)! \equiv (p-1)! \pmod p$. ✅
3. Da $(p-1)! \not\equiv 0 \pmod p$, folgt $a^{p-1} \equiv 1 \pmod p$. ✅

**Bewertung:** Vollständig korrekt. ✅

---

## Behobener Bug (BUG-B2-09)

**Bemerking nach Satz 4.1 (Wilson-Primzahlen):**

**Ursprünglicher Text (FALSCH/IRREFÜHREND):**
> A prime $p$ is a Wilson prime if and only if the even stronger congruence modulo $p$ on the Wilson quotient holds.

**Problem:** Diese Aussage ist vage und impliziert fälschlicherweise, Wilson-Primzahlen hingen mit dem vollen harmonischen Sum $H \equiv 0 \pmod p$ zusammen (was immer gilt, also trivial ist).

**Korrigierter Text (KORREKT):**
Explizite Definition der Wilson-Primzahl via $W_p = ((p-1)!+1)/p$ und Verweis auf die korrekte Formel aus Paper 10:
$$p \text{ ist Wilson-Primzahl} \iff \mu_p + S_p \equiv 0 \pmod p$$
wobei $S_p = \sum_{j=1}^{(p-1)/2} j^{-1}$ die *halbe* harmonische Summe ist (nicht die volle $H$, die stets $\equiv 0$ ist).

---

## Bugs

| ID | Status | Beschreibung |
|----|--------|-------------|
| BUG-B2-09 | ✅ Behoben | Bemerkung zu Wilson-Primzahlen: vage/irreführende Aussage durch korrekte $\mu_p+S_p$-Formel ersetzt |

---

## Fazit

Paper 12 präsentiert vier klassische Anwendungen des Wilson-Satzes mit vollständigen, korrekten Beweisen. Der einzige Fehler war eine irreführende Bemerkung zur Verbindung zwischen Wilson-Primzahlen und harmonischen Summen, die nun durch die korrekte Formel aus Paper 10 ersetzt wurde. Die Beweise für quadratische Reste und Fermats kleinen Satz sind besonders elegant.
