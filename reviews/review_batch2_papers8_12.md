# Gutachter-Bericht: Batch 2 — Wilson-Theorem-Serie (Papers 8–12)

**Datum:** 2026-03-11
**Gutachter:** Claude Sonnet 4.6
**Build:** 5
**Sprachen:** Englisch + Deutsch

---

## Überblick

Batch 2 umfasst fünf Themen zur Wilson-Theorem-Theorie, jeweils als englisches Original und deutsche Übersetzung vorhanden. Alle zehn Dateien wurden auf mathematische Korrektheit, Stimmigkeit, Vollständigkeit, logische Struktur und formale Beweis-Eigenschaften geprüft.

---

## Paper 8 — Wilson's Theorem (Grundsatz)
**Dateien:** `paper8_wilson_theorem.tex` (EN), `paper8_wilson_satz_de.tex` (DE)
**Build:** 62

### Inhalt
Vollständiger Beweis des Satzes von Wilson: $p$ prim $\Leftrightarrow$ $(p-1)! \equiv -1 \pmod{p}$.

### Mathematische Prüfung

**Struktur:**
- Lemma über Selbst-Inversen: In $(\mathbb{Z}/p\mathbb{Z})^*$ sind $\pm 1$ die einzigen Elemente mit $x^2 \equiv 1 \pmod{p}$ — korrekt, da $x^2 - 1 = (x-1)(x+1)$ in $\mathbb{F}_p$ (Integritätsbereich) genau $x = \pm 1$ als Nullstellen hat.
- Kernargument: Alle Elemente außer $\pm 1$ heben sich zu inversenPaaren auf → Produkt = $(-1) \cdot 1 = -1$ bei $p > 2$. Logisch einwandfrei.
- Rückrichtung ($n$ zusammengesetzt → $(n-1)! \not\equiv -1$): Fallunterscheidung korrekt: $n = 4$ separat behandelt ($(4-1)! = 6 \equiv 2 \pmod{4}$), allgemeines $n = a \cdot b$ mit $a \neq b$: sowohl $a$ als auch $b$ im Produkt vorhanden → $n \mid (n-1)!$. Für $n = p^2$: $p$ und $2p$ beide in $1,\ldots,p^2-1$ → $p^2 \mid (p^2-1)!$.
- Erweiterung auf Primzahlpotenzen: Korrekte Referenz auf Paper 9.

**Vollständigkeit:** ✅ Alle Fälle abgedeckt.
**Korrektheit:** ✅ Keine Fehler gefunden.
**Stimmigkeit EN ↔ DE:** ✅ Deutsche Version ist korrekte Übersetzung.

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 9 — Wilson für Primzahlpotenzen
**Dateien:** `paper9_wilson_prime_powers.tex` (EN), `paper9_wilson_primzahlpotenzen_de.tex` (DE)
**Build:** 63

### Inhalt
Beweist Wilson-Analoga für $(\mathbb{Z}/p^k\mathbb{Z})^*$: Produkt der Einheiten $\equiv -1 \pmod{p^k}$ für ungerades $p$, und $+1 \pmod{2^k}$ für $k \geq 3$.

### Mathematische Prüfung

**Zykliziät von $(\mathbb{Z}/p^k\mathbb{Z})^*$ für ungerades $p$:**
- Beweis per Hensel-Lifting: primitives Element $\bmod p$ → Lift zu primitiver Wurzel $\bmod p^k$. Korrekte Verwendung von Hensel's Lemma.
- Expliziter Konstruktionsschritt: Falls $g$ primitiv $\bmod p$, dann $g$ oder $g+p$ ist primitiv $\bmod p^2$. Korrekt.

**Struktur von $(\mathbb{Z}/2^k\mathbb{Z})^*$ für $k \geq 3$:**
- Korrekte Beschreibung: $(\mathbb{Z}/2^k\mathbb{Z})^* \cong \mathbb{Z}/2 \times \mathbb{Z}/2^{k-2}$, erzeugt von $-1$ und $5$.
- Beweis dass $\text{ord}(5) = 2^{k-2}$: per $5^{2^{k-2}} \equiv 1 \pmod{2^k}$ und $5^{2^{k-3}} \not\equiv 1 \pmod{2^k}$. Korrekt.

**Produktformel:**
- Zyklische Gruppe $\mathbb{Z}/n\mathbb{Z}$: Produkt aller Elemente = einzige Involution = Generator der Ordnung 2, d.h. $-1$. Für $p^k$ (ungerade): Produkt $\equiv -1 \pmod{p^k}$. ✅
- Für $2^k$ ($k \geq 3$): Einzige Involution in $\mathbb{Z}/2 \times \mathbb{Z}/2^{k-2}$ ist $(-1, 1)$, also Produkt $\equiv -1 \pmod{2^k}$? Hier ist Vorsicht geboten: das Paper behandelt das Produkt *aller Einheiten*, nicht nur der Gruppe als zyklische Gruppe. Die Analyse basiert auf der allgemeinen Wilson-Formel für abelsche Gruppen (Paper 11). **Konsistenz überprüft: Korrekt**, da das Produkt aller Einheiten in einer abelschen Gruppe das Produkt aller Involutionen ist (Fälle ohne eindeutige Involution ergeben $e$).

**Vollständigkeit:** ✅
**Korrektheit:** ✅
**Stimmigkeit EN ↔ DE:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 10 — Wilson-Quotient und Wilson-Primzahlen
**Dateien:** `paper10_wilson_quotient.tex` (EN), `paper10_wilson_quotient_de.tex` (DE)
**Build:** 64

### Inhalt
Definiert Wilson-Quotienten $W_p = ((p-1)!+1)/p$, beweist harmonische Kongruenz $W_p \equiv \mu_p + S_p \pmod{p}$, Bernoulli-Verbindung, Heuristik für Wilson-Primzahlen.

### Mathematische Prüfung

**Halbfaktorialformel:**
- $\mu_p = (-1)^m (m!)^2 / (-1)$ mit $m = (p-1)/2$: Korrekte Ableitung über $j \cdot (p-j) \equiv -j^2 \pmod{p}$.
- Harmonische Summe $S_p = \sum_{j=1}^{m} j^{-1} \pmod{p}$: Korrekt.
- Hauptkongruenz $W_p \equiv \mu_p + S_p \pmod{p}$: Beweis korrekt, Zwischenschritte vollständig.

**Bernoulli-Verbindung:**
- Behauptung: $W_p \equiv -B_{p-1} \pmod{p}$ (Wolstenholme-Wilson-Bernoulli-Verbindung)
- ⚠️ **SCHWÄCHE:** Der Beweis skizziert die Verbindung über $\sum_{j=1}^{p-1} j^{-1} \equiv -B_{p-1} \pmod{p}$ (via Kummer'sche Kongruenzen) ohne diese Kongruenz explizit herzuleiten. Die Gleichung $\sum_{j=1}^{p-1} j^{-1} \equiv B_{p-2} \cdot (p-1) \pmod{p}$ und der Schritt zu $B_{p-1}$ ist als "bekannte Kongruenz" zitiert. Das ist mathematisch nicht falsch, aber für ein vollständiges Paper sollte diese Verbindung referenziert oder bewiesen werden. **Empfehlung:** Explizite Quelle angeben (z.B. Granville's survey "Arithmetic properties of binomial coefficients").

**Heuristik Wilson-Primzahlen:**
- Wahrscheinlichkeitsargument: $\Pr[p \mid W_p] \approx 1/p$, erwartete Anzahl bis $x$ ≈ $\ln \ln x$. Korrekt als Heuristik.
- Bekannte Wilson-Primzahlen: 5, 13, 563 — korrekt und vollständig (Stand Wissensdatum).

**LaTeX-Bug (behoben in dieser Revision):**
- **BUG-P10-DE-LATEX:** `paper10_wilson_quotient_de.tex` verwendete `\begin{beispiel}...\end{beispiel}` ohne Definition der `beispiel`-Umgebung im Präambel. **BEHOBEN:** `\newtheorem{beispiel}[theorem]{Beispiel}` wurde ergänzt.

**Stimmigkeit EN ↔ DE:** ✅ (nach Bug-Fix)

### Urteil
**DRUCKREIF** (nach LaTeX-Fix) — Bernoulli-Beweisskizze ist inhaltlich akzeptabel für Survey-Paper, sollte aber Referenz ergänzen.

---

## Paper 11 — Wilson für endliche abelsche Gruppen
**Dateien:** `paper11_wilson_abelian.tex` (EN), `paper11_wilson_abelsch_de.tex` (DE)
**Build:** 65

### Inhalt
Verallgemeinertes Wilson-Theorem: In einer endlichen abelschen Gruppe ist das Produkt aller Elemente gleich dem Produkt aller Involutionen (Elemente der Ordnung 2); insbesondere $-1$ falls es genau eine Involution gibt, sonst $e$.

### Mathematische Prüfung

**Paarungsargument:**
- Jedes Element $g \neq g^{-1}$ (d.h. $g^2 \neq e$) paart sich mit seinem Inversen → Produkt dieser Paare = $e$. Korrekt.
- Verbleibendes Produkt = Produkt aller Involutionen (inkl. $e$). Korrekt.

**Gauss-Verallgemeinerung:**
- Theorem: $\prod_{(a,n)=1, 1 \leq a \leq n} a \equiv -1 \pmod{n}$ gdw. $n \in \{1, 2, 4, p^k, 2p^k\}$ ($p$ ungerade Primzahl, $k \geq 1$).
- Beweis via CRT-Struktursatz: $(\mathbb{Z}/n\mathbb{Z})^* \cong \prod_i (\mathbb{Z}/p_i^{k_i}\mathbb{Z})^*$.
- Genau dann eindeutige Involution, wenn das Produkt $-1$ ist. Korrekte Fallunterscheidung.
- Vollständiger Beweis aller Fälle: $n = 1, 2$: trivial; $n = 4$: Gruppe $\mathbb{Z}/2$ → einzige Involution. $n = p^k$ (ungerade): zyklische Gruppe → eindeutige Involution. $n = 2p^k$: Produkt $\cong \mathbb{Z}/2 \times (\mathbb{Z}/p^k\mathbb{Z})^*$ mit eindeutiger Involution. Übrige: Produkt enthält $\mathbb{Z}/2 \times \mathbb{Z}/2$ → kein eindeutiges $-1$.

**Vollständigkeit:** ✅ Alle Gruppen-Fälle abgedeckt.
**Korrektheit:** ✅
**Stimmigkeit EN ↔ DE:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 12 — Anwendungen des Wilson-Theorems
**Dateien:** `paper12_wilson_applications.tex` (EN), `paper12_wilson_anwendungen_de.tex` (DE)
**Build:** 66

### Inhalt
Vier Anwendungen: (1) Primzahltest, (2) Legendresymbol $\left(\frac{-1}{p}\right)$, (3) Wolstenholme harmonisches Lemma, (4) Fermats kleiner Satz via Wilson.

### Mathematische Prüfung

**Anwendung 1 — Primzahltest:**
- $p$ prim ↔ $(p-1)! \equiv -1 \pmod{p}$. Korrekt als theoretisches Kriterium.
- Korrekte Anmerkung zur Impraktikabilität ($O(p)$ Multiplikationen).

**Anwendung 2 — $\left(\frac{-1}{p}\right) = (-1)^{(p-1)/2}$:**
- Beweis: $(-1)^{(p-1)/2} \equiv \prod_{j=1}^{(p-1)/2} (-j) / \prod_{j=1}^{(p-1)/2} j = (p-1)! / ((p-1)/2)!^2 \cdot (-1)^{(p-1)/2}$...
- Das Paper verwendet den Standardbeweis: $1 \cdot (-1) \cdot 2 \cdot (-2) \cdots m \cdot (-m) = (-1)^m (m!)^2$ mit $m = (p-1)/2$, und $(p-1)! \equiv -1$ → Quadratwurzel-Argument. Korrekt.

**Anwendung 3 — Wolstenholme harmonisches Lemma:**
- Behauptung: $\sum_{j=1}^{p-1} 1/j \equiv 0 \pmod{p}$ für $p \geq 3$ prim.
- Beweis: Paare $j$ und $p-j$ ergeben $1/j + 1/(p-j) = p/(j(p-j)) \equiv 0 \pmod{p}$. Korrekt.
- **Hinweis:** Wolstenholme's starkes Resultat ($\equiv 0 \pmod{p^2}$ für $p \geq 5$) ist hier nicht behauptet, nur die schwächere Version — das ist klar ausgewiesen. ✅

**Anwendung 4 — Fermats kleiner Satz via Wilson:**
- Ableitung über Gruppentheorie: $|G| = p-1$, Lagrange → $a^{p-1} \equiv 1$. Korrekt.
- Alternativbeweis: Direkte Multiplikation $a \cdot (a \cdot 1) \cdot (a \cdot 2) \cdots = a^{p-1} (p-1)!$. Korrekt.

**Vollständigkeit:** ✅
**Korrektheit:** ✅
**Stimmigkeit EN ↔ DE:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Zusammenfassung Batch 2

| Paper | Titel | EN | DE | Urteil |
|---|---|---|---|---|
| 8 | Wilson-Grundsatz | ✅ | ✅ | DRUCKREIF |
| 9 | Primzahlpotenzen | ✅ | ✅ | DRUCKREIF |
| 10 | Wilson-Quotient | ✅ | ✅ (nach Fix) | DRUCKREIF |
| 11 | Abelsche Gruppen | ✅ | ✅ | DRUCKREIF |
| 12 | Anwendungen | ✅ | ✅ | DRUCKREIF |

### Behobene Fehler in dieser Revision
- **BUG-P10-DE-LATEX:** `beispiel`-Umgebung undefiniert → `\newtheorem{beispiel}` ergänzt.

### Offene Empfehlungen (keine Blockierer)
- Paper 10: Explizite Referenz für $\sum j^{-1} \equiv -B_{p-1} \pmod{p}$ ergänzen (z.B. Granville's survey).
