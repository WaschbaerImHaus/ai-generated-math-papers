# Gutachten: Paper 8 — Der Satz von Wilson und seine Verallgemeinerungen

**Papers:** `papers/batch2/paper8_wilson_theorem.tex` (EN) · `papers/batch2/paper8_wilson_satz_de.tex` (DE)
**Gutachter:** Michael Fuhrmann (Self-Review)
**Datum Erstgutachten:** 2026-03-11
**Datum Revisionsgutachten:** 2026-03-11
**Gesamturteil:** ✅ **Annahme ohne Auflagen** — alle 7 Bugs behoben (Build 63)

---

## Revisionsgutachten (Build 63)

### BUG-B2-01 (EN) — ✅ BEHOBEN
Abstract: „Wolstenholme primes" → „Wilson primes" korrigiert.

### BUG-B2-02 (EN+DE) — ✅ BEHOBEN
Theorem 4.1: Falsche Aussage für $p=2$, $k\geq3$ behoben.
- Vorher: Produkt $\equiv -1 \pmod{2^k}$ für $k\geq3$ (FALSCH)
- Nachher: Produkt $\equiv +1 \pmod{2^k}$ für $k\geq3$, da $(\mathbb{Z}/2^k\mathbb{Z})^*$ nicht zyklisch ist (drei Elemente der Ordnung 2)
- Korrekte Aussage: $-1$ nur für $k=1$ und $k=2$

### BUG-B2-03 (EN+DE) — ✅ BEHOBEN
Section 3 Beweis unvollständig für $n = a^2$ (Quadratzahlen ≥ 9):
- Vorher: nur Fall $a < b$ behandelt
- Nachher: expliziter Fall 2 ($n=a^2$, $a\geq3$): $a$ und $2a$ liegen in $\{1,\ldots,n-1\}$, $a\cdot2a=2n$, also $n\mid(n-1)!$

### BUG-B2-04 (EN) — ✅ BEHOBEN
Korollar 3.2: „or $n=1$" entfernt (war bei Voraussetzung $n\geq2$ sinnlos).

### BUG-B2-05 (EN) — ✅ BEHOBEN
General Wilson Beweis: `\Fp[2]` ($=\mathbb{F}_p[2]$, Polynomring) ersetzt durch korrekte Erklärung via Koordinatensumme.

### BUG-B2-06 (EN+DE) — ✅ BEHOBEN
Section 4 war nur Beweisskizze. Ersetzt durch vollständigen Beweis mit eingebettetem Lemma über das Produkt aller Elemente einer zyklischen Gruppe.

### BUG-B2-07 (DE) — ✅ BEHOBEN
Tippfehler: „Beweissskizze" → durch vollständigen Beweis ersetzt (damit obsolet).

---

## Mathematische Gesamtprüfung

| Komponente | Status | Anmerkung |
|---|---|---|
| Wilson Klassisch (if-Richtung) | ✅ | Paarungsargument korrekt |
| Wilson Klassisch (only-if-Richtung) | ✅ | Beide Fälle ($a<b$ und $a=b$) korrekt |
| Theorem 3.1 (composite $(n-1)!\equiv0$) | ✅ | Nach Fix: Fall $n=a^2$ explizit |
| Theorem 4.1 (Wilson für $p^k$) | ✅ | Nach Fix: vollständiger Beweis, $p=2$ korrekt |
| Wilson-Quotient Definition | ✅ | Korrekt |
| Wilson-Primzahlen 5, 13, 563 | ✅ | Numerisch verifiziert |
| Allgemeiner Wilson-Satz | ✅ | Drei Fälle korrekt; $|S|\geq4$-Argument korrigiert |
| Literatur | ✅ | Hardy-Wright, Ireland-Rosen, Crandall-Dilcher-Pomerance |

**Empfehlung: Druckreif.**
