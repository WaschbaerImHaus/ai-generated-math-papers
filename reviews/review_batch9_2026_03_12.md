# Review: Batch 9 — Papers 37–38 (EN + DE)
**Datum:** 2026-03-12
**Reviewer:** Mathematisches Selbst-Audit
**Papers:** paper37 (Algebraische Zahlentheorie), paper38 (Iwasawa-Theorie)

---

## Ergebnis-Übersicht

| Paper | Bugs HOCH | Bugs MITTEL | Bugs NIEDRIG | Status |
|-------|-----------|-------------|--------------|--------|
| paper37 EN | 0 | 1 | 0 | FIX NÖTIG |
| paper37 DE | 0 | 1 | 0 | FIX NÖTIG |
| paper38 EN | 0 | 1 | 0 | FIX NÖTIG |
| paper38 DE | 0 | 1 | 0 | FIX NÖTIG |

---

## Paper 37 — Algebraische Zahlentheorie (EN + DE)

### Mathematisch korrekte Teile

- **Def. Ganzheitsring / Ring of integers**: Korrekt. $\mathcal{O}_K$ als freier $\Z$-Modul vom Rang $n = [K:\Q]$. ✓
- **Quadratische Körper Integralbasen**: $d \equiv 1 \pmod 4 \Rightarrow \mathcal{O}_K = \Z[\frac{1+\sqrt{d}}{2}]$; $d \equiv 2,3 \pmod 4 \Rightarrow \mathcal{O}_K = \Z[\sqrt{d}]$. ✓
- **Diskriminante**: $\disc(\Q(\sqrt{d})) = d$ für $d \equiv 1 \pmod 4$, $4d$ sonst. ✓
- **Dedekind-Definitionsbedingungen**: Noethersch, ganzabgeschlossen, alle Nicht-Null-Primideale maximal. ✓
- **efg-Identität**: $\sum e_i f_i = n$. ✓
- **Minkowski-Schranke**: $M_K = \frac{n!}{n^n}\left(\frac{4}{\pi}\right)^{r_2}\sqrt{|\disc(K)|}$. Numerisch geprüft für $\Q(\sqrt{-5})$: $\frac{2}{\pi}\sqrt{20} \approx 2.85$ ✓
- **Dirichlet-Einheitensatz**: $\mathcal{O}_K^\times \cong \mu(K) \times \Z^{r_1+r_2-1}$. Beweis (Logarithmus-Abbildung, Kern = $\mu(K)$, Bild = volles Gitter) korrekt. ✓
- **Fundamentaleinheit $\Q(\sqrt{2})$**: $\varepsilon = 1+\sqrt{2}$, $N(\varepsilon) = -1$. ✓
- **Stark-Heegner**: Liste $d \in \{-1,-2,-3,-7,-11,-19,-43,-67,-163\}$ korrekt. ✓
- **$\Q(\sqrt{-23})$ Klassengruppe**: $h = 3$, $\mathrm{Cl}(K) \cong \Z/3\Z$. Normform $a^2+ab+6b^2$ korrekt. ✓
- **$\Q(\sqrt{5})$**: $\phi = \frac{1+\sqrt{5}}{2}$, $N(\phi) = -1$, $h = 1$ (Minkowski-Schranke $< 2$). ✓
- **Klassenzahlformel**: $\lim_{s\to 1}(s-1)\zeta_K(s) = \frac{2^{r_1}(2\pi)^{r_2} h_K R_K}{w_K \sqrt{|\disc(K)|}}$. ✓

---

### BUG-B9-P37-EN-001 — MITTEL — Theorem 6.1: Widersprüchliche Hypothese

**Datei:** `papers/batch9/paper37_algebraic_number_theory_en.tex`
**Zeilen:** Theorem (Splitting via Legendre symbol), Section 6
**Problem:**
Das Theorem beginnt mit "An odd prime $p$ **not dividing $D$**" — schließt also Ramifizierung aus.
Aber Fall 2 in der Fallunterscheidung ist "$p^2$, if $p \mid D$ (ramified)".
Das ist ein direkter Widerspruch: die Hypothese schließt den Fall aus, den Fall 2 beschreibt.

**Mathematisch korrekte Formulierung:**
Der Standardsatz gilt für ALLE ungeraden Primzahlen $p$:
- $p \nmid D$ und $(D/p) = +1$: zerfallend
- $p \mid D$: verzweigt
- $p \nmid D$ und $(D/p) = -1$: träge

**Fix:** "An odd prime $p$ not dividing $D$" → "An odd prime $p$" (Hypothese ohne $p \nmid D$).

---

### BUG-B9-P37-DE-001 — MITTEL — Satz 6.1: Gleicher Fehler in DE

**Datei:** `papers/batch9/paper37_algebraische_zahlentheorie_de.tex`
**Problem:** Identischer Widerspruch: "Eine ungerade Primzahl $p \nmid D$ erfüllt: ... falls $p \mid D$ (verzweigt)".
**Fix:** "$p \nmid D$" → kein Ausschluss, Hypothese lautet "Eine ungerade Primzahl $p$".

---

## Paper 38 — Iwasawa-Theorie (EN + DE)

### Mathematisch korrekte Teile

- **Iwasawa-Wachstumsformel**: $v_p(|A_n|) = \mu p^n + \lambda n + \nu$. ✓
- **Iwasawa-Algebra**: $\Lambda = \Zp\llbracket T \rrbracket$, vollständiger lokaler Ring, regulär, UFD, Noethersch. ✓
- **Weierstraß-Vorbereitungssatz**: $f = p^\mu \cdot u \cdot P$ mit distinguished Polynom $P$. ✓
- **Struktur-Satz für $\Lambda$-Moduln**: Pseudo-Isomorphie + charakteristisches Ideal. ✓
- **Hauptvermutung (Mazur-Wiles 1984 / Wiles 1990)**: Korrekt als bewiesener Satz deklariert. ✓
- **Ferrero-Washington**: $\mu = 0$ für abelsche Erweiterungen von $\Q$. ✓
- **Kato/Skinner-Urban**: Einseitige Teilbarkeit und $L_p(E,1) = 0 \Rightarrow \mathrm{rank}(E(\Q)) \ge 1$. Als bedingte Aussage korrekt formuliert. ✓
- **Beispiel $p=37$**: $37 \mid B_{32}$, kleinste irreguläre Primzahl, $\mu = 0$, $\lambda = 2$. ✓

---

### BUG-B9-P38-EN-001 — MITTEL — Theorem 5.1: Interpolationsbedingung zu restriktiv

**Datei:** `papers/batch9/paper38_iwasawa_theory_en.tex`
**Zeilen:** Theorem 5.1 (Kubota-Leopoldt)
**Problem:**
Die Formel wird nur für "$n \equiv 1 \pmod{p-1}$" behauptet.
Laut dem Standard-Satz (Washington: *Introduction to Cyclotomic Fields*, Theorem 5.11) gilt:
$$L_p(1-n, \chi) = (1 - \chi(p)\omega^{-n}(p)p^{n-1}) L(1-n, \chi\omega^{-n})$$
für **alle positiven ganzen Zahlen $n \geq 1$** — keine Kongruenzbedingung.
Die Einschränkung auf $n \equiv 1 \pmod{p-1}$ macht den Satz schwächer als nötig und ist für die $p$-adische Eindeutigkeit nicht ausreichend (man braucht Interpolation an unendlich vielen Stellen, die dicht liegen).

**Fix:** "$n \equiv 1 \pmod{p-1}$" → "all positive integers $n \geq 1$".

---

### BUG-B9-P38-DE-001 — MITTEL — Satz 5.1: Gleicher Fehler in DE

**Datei:** `papers/batch9/paper38_iwasawa_theorie_de.tex`
**Problem:** Identische Einschränkung "für alle positiven ganzen Zahlen $n \equiv 1 \pmod{p-1}$".
**Fix:** Kongruenzbedingung entfernen → "für alle positiven ganzen Zahlen $n \geq 1$".

---

## Zusammenfassung

| Bug-ID | Datei | Schwere | Beschreibung |
|--------|-------|---------|--------------|
| BUG-B9-P37-EN-001 | paper37_en | MITTEL | Theorem 6.1: Hypothese "$p \nmid D$" widerspricht Fall "$p \mid D$" |
| BUG-B9-P37-DE-001 | paper37_de | MITTEL | Gleicher Fehler in DE |
| BUG-B9-P38-EN-001 | paper38_en | MITTEL | Kubota-Leopoldt: Interpolation zu restriktiv ($n \equiv 1 \pmod{p-1}$) |
| BUG-B9-P38-DE-001 | paper38_de | MITTEL | Gleicher Fehler in DE |

**Keine kritischen oder hohen Bugs.** Alle 4 MITTEL-Bugs sind klar behebbar.
