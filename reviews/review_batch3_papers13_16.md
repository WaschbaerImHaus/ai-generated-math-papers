# Gutachter-Bericht: Batch 3 — Sieblehre-Serie (Papers 13–16)

**Datum:** 2026-03-11
**Gutachter:** Claude Sonnet 4.6
**Build:** 5
**Sprachen:** Englisch + Deutsch

---

## Überblick

Batch 3 umfasst vier Themen der analytischen Sieblehre, jeweils als englisches Original und deutsche Übersetzung vorhanden. Geprüft auf mathematische Korrektheit, Stimmigkeit, Vollständigkeit und Beweiseigenschaften.

---

## Paper 13 — Bruns Satz
**Dateien:** `paper13_bruns_theorem.tex` (EN), `paper13_brunnscher_satz_de.tex` (DE)
**Build:** 74

### Inhalt
Vollständiger Beweis von Bruns Satz: $\sum_{p, p+2 \text{ Primzwillinge}} (1/p + 1/(p+2)) < \infty$ (Brun-Konstante $B_2 \approx 1.902$).

### Mathematische Prüfung

**Kombinatorisches Sieb (Abschneidungsmethode):**
- Inklusionsexklusion mit Abschneidung auf gerade Stufe: $S(A, z) \leq \sum_{d \mid P(z), \omega(d) \leq 2k} \mu(d) |A_d|$ (Bonferroni-Ungleichung). Korrekt.
- Abschätzung der Restterme: $|A_d| = X/d + O(1)$ → Hauptterm und Fehlerterm getrennt behandelt. Korrekt.

**Obere Schranke für $\pi_2(x)$:**
- $\pi_2(x) = O\!\left(\frac{x (\log \log x)^2}{(\log x)^2}\right)$: Korrekte Aussage.
- Beweis: Hauptterm-Abschätzung $\prod_{p \leq z}(1 - 2/p) \approx C/(\log z)^2$ via Mertens' drittes Theorem. ✅

**Konvergenz der Zwillingsprimzahlreziproken:**
- Dyadische Partialsummenabschätzung: Summation über Blöcke $[2^k, 2^{k+1})$, jeder Block trägt $O(1/k^2)$ bei → Gesamtreihe konvergent. Korrekt.

**Numerische Tabelle:**
- Bekannte Schranken und numerische Werte für $B_2$ korrekt angegeben.

**Vollständigkeit:** ✅
**Korrektheit:** ✅
**Stimmigkeit EN ↔ DE:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 14 — Selbergs Sieb
**Dateien:** `paper14_selberg_sieve.tex` (EN), `paper14_selberg_sieb_de.tex` (DE)
**Build:** 74

### Inhalt
Vollständige Darstellung von Selbergs $\lambda^2$-Methode, Diagonalisierung via Möbius-Inversion, Brun-Titchmarsh-Theorem.

### Mathematische Prüfung

**$\lambda^2$-Methode (quadratische Optimierung):**
- Ansatz: $\sum_{n \in A, (n, P(z))=1} 1 \leq \sum_{n \in A} \left(\sum_{d \mid (n, P(z))} \lambda_d\right)^2$ mit $\lambda_1 = 1$. Korrekt als obere Schranke.
- Optimierung über $\lambda_d$: Führt zu $V(z) = \sum_{d \mid P(z)} \mu^2(d) \prod_{p \mid d} \frac{g(p)}{1-g(p)}$. Explizite Ableitung per Diagonalisierung. ✅

**Möbius-Inversion:**
- Schlüsselidentität $V_h(z) = 1/V(z)$: Die Herleitung über die Diagonalform der quadratischen Form ist vollständig und korrekt dargestellt.

**Brun-Titchmarsh-Theorem:**
- Behauptung: $\pi(x; q, a) \leq \frac{2x}{\phi(q) \log(x/q)}$ für $x > q$. Korrekt.
- Beweis durch Anwendung des Selberg-Siebs auf $A = \{n \leq x : n \equiv a \pmod{q}\}$. Vollständig.

**Parität-Problem (Parity Problem):**
- Korrekte Diskussion der Selberg-Schranke: Unterschranken nicht direkt erreichbar (Parity obstruction). Referenz auf Chen's Theorem als Überwindung durch Buchstab-Umschalten. ✅

**Vollständigkeit:** ✅
**Korrektheit:** ✅
**Stimmigkeit EN ↔ DE:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 15 — Die große Sieb-Ungleichung
**Dateien:** `paper15_large_sieve.tex` (EN), `paper15_grosses_sieb_de.tex` (DE)
**Build:** 56 → **74 korrigiert** (Build-Inkonsistenz behoben)

### Inhalt
Die Montgomery-Vaughan-Form der großen Sieb-Ungleichung: $\sum_{q \leq Q} \sum_{\substack{a=1\\(a,q)=1}}^{q} \left|\sum_{n=M+1}^{M+N} a_n e(an/q)\right|^2 \leq (N + Q^2) \sum_{n=M+1}^{M+N} |a_n|^2$.

### Mathematische Prüfung

**Fourier-analytischer Kern:**
- Beweis basiert auf dem Fejér-Kern-Ansatz: Nutzung der Orthogonalitätseigenschaften der Farey-Folge. Korrekt im Ansatz.
- Die Abschätzung $\delta$-well-spaced sets: Punkte $\alpha_r$ mit $\|\alpha_r - \alpha_s\| \geq \delta$ für $r \neq s$ → Hauptterm-Dominanz. ✅

**Beweisskizze (korrekt als solche gekennzeichnet):**
- Das Paper markiert den vollständigen Beweis ausdrücklich als "Proof Sketch" und verweist auf Montgomery (1971). Dies ist methodisch korrekt für ein Survey-Paper.

**Anwendung: Bombieri-Vinogradov:**
- Korrekte Aussage: $\sum_{q \leq x^{1/2}/(\log x)^A} \max_{(a,q)=1} \left|\psi(x;q,a) - \frac{x}{\phi(q)}\right| = O\!\left(\frac{x}{(\log x)^B}\right)$.
- Beweis-Idee korrekt skizziert (GRH im Durchschnitt).

**Anwendung: Nullstellenfreie Region:**
- Korrekte Verknüpfung mit dem Null-Dichte-Theorem (zero-density estimate).

**Build-Fehler (behoben):**
- **BUG-BUILD-INCONSISTENCY:** Build 56 statt korrektem Build 74. **BEHOBEN.**

**Vollständigkeit:** ✅ (als Survey-Paper, Beweisskizze klar ausgewiesen)
**Korrektheit:** ✅
**Stimmigkeit EN ↔ DE:** ✅

### Urteil
**DRUCKREIF** — keine Beanstandungen.

---

## Paper 16 — Chens Satz
**Dateien:** `paper16_chen_theorem.tex` (EN), `paper16_chen_satz_de.tex` (DE)
**Build:** 56 → **74 korrigiert** (Build-Inkonsistenz behoben)

### Inhalt
Chen's Theorem (1966/1973): Jede hinreichend große gerade Zahl $N = p + q$ mit $p$ prim und $q \in \mathcal{P}_2$ (Primzahl oder Semiprimes).

### Mathematische Prüfung

**Beweis-Strategie:**
- Zwei-Schritt-Ansatz: Obere Schranke für $\mathcal{P}_2$-Terme + Buchstab-Umschaltprinzip für Untere Schranke. Korrekt.

**Ingredient 1 — Selberg Oberschranke:**
- Anwendung des Selberg-Siebs auf $A = \{N - p : p \leq N, p \text{ prim}\}$. Korrekte Anwendung.

**Ingredient 2 — Rosser-Iwaniec Unterschranke:**
- Nutzung des linearen Buchstab-Siebs mit Iwaniec-Verbesserungen. Korrekt.
- Der Faktor $\frac{1}{2}$ aus dem Buchstab-Umschaltprinzip ist korrekt erklärt: Das Switching-Prinzip reduziert den Koeffizienten der unteren Schranke von $\sim 1$ auf $\sim \frac{1}{2}$, was für Chens Beweis ausreicht.

**Parität-Problem:**
- Korrekter Hinweis, dass Chen's Beweis das Parity-Problem umgeht (nicht löst) durch spezifische Gewichtsmethode.

**Behobene Fehler in dieser Revision:**
- **BUG-P16-DE-MISSING-REF:** `paper16_chen_satz_de.tex` fehlte `\bibitem{Iwaniec1980}` (H. Iwaniec, "Rosser's sieve", Acta Arith. 36, 1980). **BEHOBEN.**
- **BUG-BUILD-INCONSISTENCY:** Build 56 → 74. **BEHOBEN.**

**Vollständigkeit:** ✅
**Korrektheit:** ✅
**Stimmigkeit EN ↔ DE:** ✅ (nach Bug-Fix)

### Urteil
**DRUCKREIF** (nach Fixes) — keine weiteren Beanstandungen.

---

## Zusammenfassung Batch 3

| Paper | Titel | EN | DE | Urteil |
|---|---|---|---|---|
| 13 | Bruns Satz | ✅ | ✅ | DRUCKREIF |
| 14 | Selbergs Sieb | ✅ | ✅ | DRUCKREIF |
| 15 | Große Sieb-Ungleichung | ✅ | ✅ (nach Fix) | DRUCKREIF |
| 16 | Chens Satz | ✅ | ✅ (nach Fixes) | DRUCKREIF |

### Behobene Fehler in dieser Revision
- **BUG-BUILD-INCONSISTENCY:** Build 56 in Papers 15 (EN+DE) und 16 (EN+DE) → korrigiert auf Build 74.
- **BUG-P16-DE-MISSING-REF:** Fehlende Bibitem `Iwaniec1980` in DE-Version von Paper 16 ergänzt.
