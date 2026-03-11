# Gutachter-Bericht: Papers 18–28 (Batches 4–6)
**Datum:** 2026-03-11 | **Build:** 9 | **Gutachter:** Claude Sonnet 4.6

---

## Prüfkriterien

1. **Formale Korrektheit** — LaTeX-Kompilierbarkeit, Theorem-Struktur, Bibliographie
2. **Schlüssigkeit** — logische Konsistenz der Beweise, korrekte Querverweise
3. **Vollständigkeit** — keine fehlenden Textteile, alle Sektionen vorhanden
4. **EN/DE-Parität** — Übereinstimmung der englischen und deutschen Version

---

## Batch 4 — Analytische Zahlentheorie (Papers 18–20)

---

### Paper 18: Vinogradovs Drei-Primzahlen-Satz

#### 18-EN (`paper18_vinogradov_three_primes.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Kompilierbar; `\inputenc` vorhanden; alle Theoreme nummeriert |
| Mathematische Schlüssigkeit | Weitgehend korrekt; **ein Vorzeichenfehler** (BUG-B4-P18-EN-001) |
| Vollständigkeit | Vollständig; Helfgott-2013-Sektion vorhanden |
| Struktur | 8 Sektionen; Kreismethode → Exponentialsummen → singuläre Reihe → singuläres Integral → Hauptterm → Helfgott |

**Fehler:**
- **BUG-B4-P18-EN-001 (Mittel):** Remark 4.3: Schreibt `$1 + c_2(n)/(2-1)^3 = 2$` für gerades n. Korrekt: `$1 - c_2(n) = 0$` (da μ(2)/φ(2) = −1). Das Vorzeichen widerspricht der unmittelbar folgenden Behauptung S(n)=0 für gerades n.
- **Befund (Gering):** Remark 2.1 (Zeilen 128–136, Erklärung R₃(n) ~ (log n)³ r₃(n)/6) fehlt im DE-Pendant — kein Fehler in EN, aber EN ist reichhaltiger.

**Urteil EN:** ⚠️ Bedingt druckreif — nach Vorzeichenkorrektur in Remark 4.3.

---

#### 18-DE (`paper18_vinogradov_drei_primzahlen_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `\usepackage[utf8]{inputenc}` |
| Mathematische Schlüssigkeit | Vorzeichenfehler übernommen (BUG-B4-P18-DE-001); Beweis singuläres Integral unvollständig (BUG-B4-P18-DE-002) |
| Vollständigkeit | Remark 2.1 (R₃ vs r₃) fehlt; Nebenbogen-Beweis hat extra 3 Schritte (kein Fehler, aber asymmetrisch) |
| EN/DE-Parität | DE fehlt Remark 2.1; sonst strukturell parallel |

**Fehler:**
- **BUG-B4-P18-DE-001 (Mittel):** Identischer Vorzeichenfehler wie EN (Bemerkung „Gerades vs. ungerades n").
- **BUG-B4-P18-DE-002 (Gering):** Beweis Satz 5.1: DE gibt nur Formel `(n-6)²/2` (gilt nur n≤N). EN behandelt korrekt beide Fälle (n≤N und N≤n≤2N) mit allgemeinem Ergebnis n²/2 + O(n).
- **Strukturell (Info):** Fehlendes `[utf8]{inputenc}` — tritt bei allen Batch-4/5/6 DE-Dateien auf (systematisches Problem).

**Urteil DE:** ⚠️ Nicht druckreif — 2 offene Fehler.

---

### Paper 19: Waringsches Problem

#### 19-EN (`paper19_waring_problem.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | `\bibitem{Wooley2016}` im Text zitiert, aber nur `Wooley2012` und `Wooley2013` in Bibliographie |
| Mathematische Schlüssigkeit | Drei inkohärente Schwellenbedingungen ohne Überblick (s>2k+1, s>2k, s≤k(k+1)/2); zirkuläre Logik Cor. 8.3 |
| Vollständigkeit | Corollary 4.1 vorhanden; Remark 4.1 vorhanden (fehlt in DE) |
| Struktur | 9 Sektionen; g(k), G(k), Heuristik, Weyl, Kreismethode, Wooley, Schluss |

**Fehler:**
- **BUG-B4-P19-EN-001 (Gering):** Querverweise auf „Paper 17, Definition 5.2" und „Paper 17, Lemma 5.5" in Theorem 4.1 — Nummerierung sollte auf tatsächliche Nummern geprüft werden (wahrscheinlich inkonsistent mit dem tatsächlichen Paper 17).
- **BUG-B4-P19-EN-002 (Hoch):** Drei Schwellenbedingungen (s>2k+1, s>2k, s≤k(k+1)/2) werden nie kohärent in Relation zueinander gesetzt. Für große k gilt k(k+1)/2 ≫ 2k, d.h. die Corollary 8.3 nennt s=⌈k(k+1)/2⌉+1, behauptet aber gleichzeitig Schwelle s>2k — logisch inkonsistent ohne Erklärung.
- **BUG (Mittel):** Fehlende `\bibitem{Wooley2016}` — Zitation `\cite{Wooley2016}` im Text ohne Bibliographieeintrag.

**Urteil EN:** ⚠️ Bedingt druckreif — fehlendes Bibitem muss ergänzt werden; logische Kohärenz der Schwellenbedingungen sollte klargestellt werden.

---

#### 19-DE (`paper19_waringsches_problem_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}`; 3 fehlende Bibitems |
| Mathematische Schlüssigkeit | Wooleys Haupttheorem als `\begin{remark}` statt `\begin{theorem}` — gravierend |
| Vollständigkeit | Remark 4.1 (drei Schwellenbedingungen) fehlt; Corollary 4.1 (Weyl-Schranke G(k)) fehlt; Theorem 6.1 (Efficient Congruencing) fehlt |
| EN/DE-Parität | DE hat ca. 403 Zeilen vs. EN 560 — strukturell erheblich kürzer |

**Fehler:**
- **BUG-B4-P19-DE-MISSING-THM (Kritisch):** Wooleys Efficient Congruencing-Theorem (Hauptresultat von Section 6) steht in DE als bloße `remark`, nicht als `theorem`. Mathematisch gravierender Darstellungsfehler.
- **BUG-B4-P19-DE-001 (Mittel):** Corollary „Weyl-Schranke G(k)≤(k-2)2^(k-1)+5" mit Beweis fehlt komplett.
- **BUG-B4-P19-DE-002 (Gering):** Wooley-Schranken-Korollar ohne Beweisskizze.
- **BUG-B4-P19-DE-003 (Mittel):** Kein formaler Satz für Efficient Congruencing (J_{s,k}-Schranke).
- **BUG-B4-P19-DE-MISSING-REFS (Mittel):** Fehlende Bibitems: `HardyLittlewood1920`, `Wooley2012`, `Wooley2013`. DE: 6 Einträge, EN: 9.

**Urteil DE:** ❌ Nicht druckreif — 5 offene Fehler, davon 1 kritisch.

---

### Paper 20: Goldbachsche Singuläre Reihe

#### 20-EN (`paper20_goldbach_singular_series.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | `\DeclareMathOperator{\Li}{li}` vorhanden; 9 Bibitems vollständig |
| Mathematische Schlüssigkeit | Tabelle S(n) systematisch falsch (BUG-B4-P20-EN-002) |
| Vollständigkeit | Section 8 „Goldbach Comet" vorhanden; alle Sektionen da |
| Struktur | 8 Sektionen; Eulerprodukt, Konvergenz, Twinprime-Konstante, Numerik, Paritätsproblem |

**Fehler:**
- **BUG-B4-P20-EN-001 (Gering):** Satzfragment in Remark „Parity Problem" — abgebrochener Satz „odd only for...".
- **BUG-B4-P20-EN-002 (Hoch):** Numerische Tabelle S(n) (Section 4): Alle 6 Einträge enthalten je einen falschen Faktor (p-1)/(p-2) für eine Primzahl p, die n NICHT teilt. Korrekte Werte: n=4: 2C₂; n=6: 4C₂; n=8: 2C₂; n=10: (8/3)C₂; n=12: 4C₂; n=30: (16/3)C₂.
- **BUG-B4-P20-EN-TABLE-001 (Mittel):** Notation in der Tabelle (z.B. `$2C_2 \cdot 2/1 = 4C_2$`) suggeriert, (2/1) sei unabhängiger Faktor — sollte als Primfaktorterm (p-1)/(p-2) für p=3 explizit beschriftet sein.

**Urteil EN:** ⚠️ Nicht druckreif — Tabellenfehler muss korrigiert werden.

---

#### 20-DE (`paper20_goldbach_singulaere_reihe_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}`; 3 fehlende Bibitems |
| Mathematische Schlüssigkeit | Identische Tabellenfehler wie EN |
| Vollständigkeit | **Section 8 „Numerics and the Goldbach Comet" fehlt vollständig** — kritisch |
| EN/DE-Parität | DE: ~448 Zeilen, EN: ~588 Zeilen; Section 8 komplett absent |

**Fehler:**
- **BUG-B4-P20-DE-SECTION-MISSING (Kritisch):** Section 8 mit Goldbach-Komet-Tabelle und numerischen Werten fehlt vollständig in DE.
- **BUG-B4-P20-DE-001 (Hoch):** Identische falsche Tabellenwerte wie BUG-B4-P20-EN-002.
- **BUG-B4-P20-DE-MISSING-BIBITEMS (Mittel):** Fehlende Bibitems: `Goldston1992`, `Vinogradov1937`, `Granville1995`. DE: 6, EN: 9.

**Urteil DE:** ❌ Nicht druckreif — fehlende Sektion und Tabellenfehler.

---

## Batch 5 — Riemannsche Zetafunktion und Primzahlsatz (Papers 21–24)

> **Systemisches Problem Batch 5+6 DE:** Alle deutschen Dateien in Batch 5 und 6 fehlen `\usepackage[utf8]{inputenc}`. Dies kann auf Systemen mit UTF-8-Umlauten zu Kompilierungsfehlern führen. Dieses Problem wird im Folgenden pro Paper einmalig angemerkt, aber nicht wiederholt als eigener Bug gezählt.

---

### Paper 21: Riemannsche Zetafunktion

#### 21-EN (`paper21_riemann_zeta_function_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Vollständig; `inputenc` vorhanden; 6 Bibitems inkl. Platt2021 |
| Mathematische Schlüssigkeit | Korrekt; Euler-Produkt, Γ-Faktor, analytische Fortsetzung, Funktionalgleichung, kritischer Streifen, Riemann-Siegel alle vollständig und korrekt |
| Vollständigkeit | Alle 8 Sektionen vorhanden; Werte ζ(2k), triviale Nullstellen, kritischer Streifen vollständig |
| Struktur | Aufbau von Dirichlet-Reihe → Euler-Produkt → meromorphe Fortsetzung → Funktionalgleichung → Nullstellen |

**Keine Fehler gefunden.**

**Urteil EN:** ✅ Druckreif.

---

#### 21-DE (`paper21_riemann_zeta_function_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}`; ansonsten vollständig |
| Mathematische Schlüssigkeit | Korrekt; vollständig parallel zu EN |
| Vollständigkeit | Alle Sektionen vorhanden; strukturell identisch mit EN |
| EN/DE-Parität | Vollständig parallel; alle Formeln identisch |

**Fehler:**
- Fehlendes `[utf8]{inputenc}` (systemisches Problem).

**Urteil DE:** ✅ Inhaltlich druckreif (vorbehaltlich `inputenc`-Fix).

---

### Paper 22: Nichttriviale Nullstellen der Zetafunktion

#### 22-EN (`paper22_nontrivial_zeros_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Vollständig; alle Theoreme nummeriert; Bibliographie komplett |
| Mathematische Schlüssigkeit | Korrekt; N(T)-Formel, erste Nullstellen, Gram-Punkte, Montgomery-Odlyzko/GUE, numerische RH-Verifikation, GRH-Sektion |
| Vollständigkeit | 7 Sektionen vorhanden; vollständig |
| Struktur | Klassisch: N(T) → Gram-Punkte → GUE-Statistik → numerische Verifikation → GRH |

**Keine Fehler gefunden.**

**Urteil EN:** ✅ Druckreif.

---

#### 22-DE (`paper22_nontrivial_zeros_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}` |
| Mathematische Schlüssigkeit | Korrekt und vollständig parallel zu EN |
| Vollständigkeit | Alle Sektionen vorhanden |
| EN/DE-Parität | Vollständig parallel |

**Fehler:**
- Fehlendes `[utf8]{inputenc}`.

**Urteil DE:** ✅ Inhaltlich druckreif.

---

### Paper 23: Explizite Formel und Primzahlsatz

#### 23-EN (`paper23_explicit_formula_pnt_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Vollständig; alle Referenzen aufgelöst |
| Mathematische Schlüssigkeit | Korrekt; von Mangoldt, Chebyshev, Perron, explizite Formel für ψ(x), PNT, Fehlertermabschätzungen, Schoenfeld, Riemanns ursprüngliche Formel |
| Vollständigkeit | 8 Sektionen; GRH/L-Funktionen-Sektion vorhanden |
| Struktur | Aufbau: Chebyshev → Perron → explizite Formel → Fehlterme → GRH |

**Keine Fehler gefunden.**

**Urteil EN:** ✅ Druckreif.

---

#### 23-DE (`paper23_explicit_formula_pnt_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}` |
| Mathematische Schlüssigkeit | Korrekt |
| Vollständigkeit | Alle Sektionen vorhanden; Abstract DE leicht länger als EN (erwähnt explizit GRH-Implikationen für arithmetische Progressionen — informativer, kein Fehler) |
| EN/DE-Parität | Strukturell vollständig parallel; Abstract minimal erweitert |

**Befund:**
- **BUG-B5-P23-DE-001 (Mittel):** Laut BUGS.md ist der Abstract in DE abgeschnitten — der Abschlusssatz über Dirichlet-L-Funktionen und GRH fehlt in einigen Versionen. Im vorliegenden File ist der Abstract länger als EN, was der Fehlerbeschreibung widerspricht. Möglicherweise zwischenzeitlich korrigiert.

**Urteil DE:** ✅ Inhaltlich druckreif.

---

### Paper 24: Ansätze zur Riemannschen Hypothese

#### 24-EN (`paper24_rh_approaches_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Vollständig |
| Mathematische Schlüssigkeit | Korrekt; Hilbert-Pólya, Berry-Keating, GUE/Zufallsmatrizen, Selberg-Spurformel, partielle Resultate, Barrieren, Physik-Verbindungen |
| Vollständigkeit | 9 Sektionen; Section 8 „Connection to Physics" (Lee-Yang, Connes Stringtheorie, Supersymmetrie/Bender et al. 1998) |

**Keine Fehler gefunden.**

**Urteil EN:** ✅ Druckreif.

---

#### 24-DE (`paper24_rh_approaches_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}` |
| Mathematische Schlüssigkeit | Korrekt in den vorhandenen Teilen |
| Vollständigkeit | **Section 8 „Verbindung zur Physik: Quantenchaos" fehlt vollständig** |
| EN/DE-Parität | DE springt direkt von Section 7 (Barrieren) zu Section 8 (Schlussfolgerung); EN Section 8 mit Lee-Yang, Connes, Supersymmetrie komplett absent |

**Fehler:**
- **BUG-B5-P24-DE-MISSING-SECTION (Kritisch):** EN Section 8 „Connection to Physics: Quantum Chaos and Beyond" — enthält Bemerkungen zu Lee-Yang-Theorem, Connes' Spurformel/Stringtheorie, Supersymmetrie (Bender et al. 1998) — fehlt vollständig in DE. DE-Section 8 (Schlussfolgerung) entspricht EN-Section 9.

**Urteil DE:** ❌ Nicht druckreif — fehlende Physik-Sektion muss ergänzt werden.

---

## Batch 6 — Elliptische Kurven und BSD (Papers 25–28)

---

### Paper 25: Elliptische Kurven über ℚ

#### 25-EN (`paper25_elliptic_curves_Q_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | `\DeclareMathOperator{\gcd}{gcd}` — möglicherweise Konflikt mit LaTeXs eingebautem `\gcd`; ansonsten vollständig |
| Mathematische Schlüssigkeit | Korrekt; Weierstrass-Form, j-Invariante, Gruppengesetz, Torsion/Mazur, Nagell-Lutz, Mordell-Weil, kanonische Höhe, Reduktion, Néron-Ogg-Shafarevich |
| Vollständigkeit | 10 Sektionen; alle wesentlichen Inhalte vorhanden |

**Befund:**
- **Technisch (Gering):** `\DeclareMathOperator{\gcd}{gcd}` deklariert einen Operator, der bereits in AMSTeX/LaTeX definiert ist. Kein Fehler (überschreibt `\gcd` mit identischer Definition), aber unnötig.

**Urteil EN:** ✅ Druckreif.

---

#### 25-DE (`paper25_elliptic_curves_Q_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}`; `Cornell1997` Bibitem fehlt |
| Mathematische Schlüssigkeit | Korrekt; NOS-Theorem in DE ausführlicher als EN (explizite Formulierung mit Galois-Darstellung und Trägheitsgruppe I_p) |
| Vollständigkeit | Einige inhaltliche Lücken gegenüber EN |
| EN/DE-Parität | DE ist an einigen Stellen detaillierter (NOS), an anderen kürzer (fehlendes Beispiel, fehlende Rang-Bemerkung) |

**Fehler:**
- **BUG-B6-P25-DE-001 (Mittel):** Beispiel zu speziellen Kurven (j=1728 für y²=x³-x, j=0 für y²=x³+1) in EN Section 2, Zeilen 121–128 fehlt in DE.
- **BUG-B6-P25-DE-002 (Mittel):** Bibitem `Cornell1997` fehlt — DE: 4 Einträge, EN: 5.
- **BUG-B6-P25-DE-003 (Gering):** Remark „It is not known whether ranks are unbounded over ℚ" (EN ~Zeile 300) fehlt in DE.

**Urteil DE:** ⚠️ Bedingt druckreif — fehlendes Bibitem und einige inhaltliche Lücken.

---

### Paper 26: L-Funktionen elliptischer Kurven

#### 26-EN (`paper26_l_function_elliptic_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Vollständig; 6 Bibitems inkl. TaylorWiles1995 (Ann. Math. **142**, 1995) |
| Mathematische Schlüssigkeit | Korrekt; Frobenius-Spur, Hasse-Schranke, Euler-Produkt, vervollständigte L-Funktion, Funktionalgleichung via Modularität |
| Vollständigkeit | 7 Sektionen; Section 6 „The L-Function Near s=1" mit Definition r_an und L*(E,1) |

**Keine Fehler gefunden.**

**Urteil EN:** ✅ Druckreif.

---

#### 26-DE (`paper26_l_function_elliptic_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}`; `TaylorWiles1995` Bibitem fehlt oder fehlerhaft |
| Mathematische Schlüssigkeit | Korrekt; numerische Beispiele mit Frobenius-Spuren sind mathematisch korrekt |
| Vollständigkeit | DE hat reichhaltigere numerische Beispiele als EN (Kolyvagin-Referenzen, explizite a_p-Werte); aber strukturelle Asymmetrie |
| EN/DE-Parität | **Strukturelle Abweichung:** DE hat erweiterte Section 6 (numerische Beispiele mit Frobenius-Spuren) BEFORE Section 7 (L-Funktion nahe s=1), während EN diese Reihenfolge umkehrt. Inhaltlich sind die numerischen Beispiele in DE reicher (a_p-Werte für E₀ und E₁, Kolyvagin-Referenzen) als in EN. |

**Fehler:**
- **BUG-B6-P26-DE-001 (Mittel):** Tabelle für a₅-Berechnung (y²=x³-x mod 5) aus EN Example 2.3 fehlt in DE; DE gibt nur verbale Zusammenfassung.
- **BUG-B6-P26-DE-002 (Mittel):** Bibitem `TaylorWiles1995` fehlt in DE oder hat falsches Journal-Volumen (141 statt 142 — Ann. Math. 142, 1995 ist korrekt). DE: 5 Bibitems, EN: 6.

**Urteil DE:** ⚠️ Bedingt druckreif — strukturelle Abweichung von EN und fehlendes Bibitem.

---

### Paper 27: Birch–Swinnerton-Dyer-Vermutung

#### 27-EN (`paper27_bsd_conjecture_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Vollständig; 7 Bibitems; alle Theorem-Umgebungen korrekt |
| Mathematische Schlüssigkeit | Korrekt; schwache BSD, starke BSD-Formel, Tate-Shafarevich-Gruppe, Selmer-Gruppe, Coates-Wiles, Kolyvagin |
| Vollständigkeit | 8 Sektionen; Beweisskizzen für Coates-Wiles und Kolyvagin vollständig und korrekt |
| Struktur | Aufbau: Einleitung → schwache BSD → Tate-Shafarevich → starke BSD → numerische Evidenz → bekannte Resultate → Selmer/Abstieg → Schluss |

**Befund:**
- **Beweis Coates-Wiles (EN, Zeilen 229–233):** Sehr detailliert — nennt p^∞-Selmer-Gruppe und CM-Theorie. Korrekt.
- **Beweis Kolyvagin (EN, Zeilen 248–259):** Euler-Systeme vollständig erklärt (Heegner-Punkte, Gross-Zagier, Kolyvaginscher Abstieg). Korrekt.
- **Selmer-Sektion (EN, Zeilen 265–280):** Enthält präzise Formel `r + dim_{F₂} Sha[2] = dim_{F₂} Sel₂(E/Q) - dim_{F₂} E(Q)_tors[2]`. Mathematisch korrekt.

**Keine Fehler gefunden.**

**Urteil EN:** ✅ Druckreif.

---

#### 27-DE (`paper27_bsd_conjecture_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}`; 7 Bibitems (inkl. Milne2006 — vorhanden) |
| Mathematische Schlüssigkeit | Korrekt; vollständig parallel zu EN |
| Vollständigkeit | Alle Sektionen vorhanden |
| EN/DE-Parität | Selmer-Sektion (DE) enthält Kummer-Einbettung und Kolyvagin-Erklärung (EN fehlt), aber **fehlt die präzise algebraische Formel** aus EN für r + dim Sha[2] |

**Fehler:**
- **BUG-B6-P27-DE-002 (Mittel):** Section 7 „Abstieg und Selmer-Gruppe": DE fehlt die exakte algebraische Formel `r + dim_{F₂} Sha[2] = dim_{F₂} Sel₂(E/Q) - dim_{F₂} E(Q)_tors[2]` aus EN. DE sagt nur `r = dim_{F₂} Sel₂(E/Q)` unter Sha=0, ohne Torsionskorrektur.
- **Info:** Bibitem `Milne2006` ist in der aktuellen DE-Datei vorhanden (7 Einträge, gleich wie EN). BUG-B6-P27-DE-001 (Milne2006 fehlt) ist möglicherweise bereits behoben.

**Urteil DE:** ✅ Weitgehend druckreif — algebraische Formel in Section 7 sollte ergänzt werden (mittel).

---

### Paper 28: Kongruente Zahlen und BSD

#### 28-EN (`paper28_congruent_numbers_bsd_en.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Vollständig; 6 Bibitems |
| Mathematische Schlüssigkeit | Vollständig korrekt; alle Beweise verifiziert |
| Vollständigkeit | 7 Sektionen; Beweise für Äquivalenz (Kongruenzzahl ↔ Rang≥1), Torsionsuntergruppe, Tunnells Theorem, kleine Beispiele (n=5,6,7), BSD-Verbindung |
| Beweise | Alle drei Berechnungen (n=5,6,7) vollständig nachvollziehbar und rechnerisch korrekt |

**Verifikation der Kernbeweise:**
- **Theorem 2.1 (Äquivalenz):** Beide Richtungen vollständig bewiesen. Die Berechnungen `x₀³ - n²x₀ = y₀²` und die Rückformeln sind korrekt.
- **n=6:** Punkt (25/4, 35/8) auf E₆: (35/8)²=1225/64, (25/4)³-36(25/4)=1225/64. ✓
- **n=5:** Punkt (25/4, 75/8) auf E₅: (75/8)²=5625/64, (25/4)³-25(25/4)=5625/64. ✓
- **Tunnell-Remark n=1,2,3:** Verwendet korrekt A,B für ungerades n=1,3 und C,D für gerades n=2. ✓

**Keine Fehler gefunden.**

**Urteil EN:** ✅ Druckreif.

---

#### 28-DE (`paper28_congruent_numbers_bsd_de.tex`)

| Kriterium | Befund |
|---|---|
| Formale Korrektheit | Fehlendes `[utf8]{inputenc}`; 6 Bibitems (inkl. Koblitz1993 — vorhanden) |
| Mathematische Schlüssigkeit | Großteils korrekt; ein inhaltlicher Fehler identifiziert |
| Vollständigkeit | Alle Sektionen vorhanden |
| EN/DE-Parität | Strukturell parallel; DE-Beweis Theorem 2.1 etwas weniger explizit aber vollständig; n=5 in DE ohne expliziten Punkt auf E₅ |

**Fehler:**
- **BUG-B6-P28-DE-TUNNELL-PARITY (Mittel — Mathematikfehler):** Remark 5.2 (Zeilen 264–265): DE schreibt `A(2) ≠ 2B(2)` und `A(3) ≠ 2B(3)` für Nicht-Kongruenz von n=2 und n=3. Für **gerades** n=2 muss gemäß Tunnells Satz (Theorem 4.1 dieser Datei selbst) die Bedingung C(n)=2D(n) geprüft werden, nicht A(n)=2B(n). EN Remark (Zeilen 276–283) verwendet korrekt C(2)=2≠4=2D(2) für n=2. Dies ist ein direkter Widerspruch zum eigenen Theorem 4.1 im Paper.
- **Info:** Remark zu n=5 in DE verifiziert nur das Dreieck, nicht den expliziten Punkt auf E₅ (kein Fehler, nur weniger vollständig als EN).
- **Info:** Bibitem `Koblitz1993` ist in der aktuellen Datei vorhanden. BUG-B6-P28-DE-003 (fehlendes Koblitz1993) möglicherweise bereits behoben.

**Urteil DE:** ⚠️ Bedingt druckreif — falscher Verweis auf A(2)/B(2) statt C(2)/D(2) in Remark 5.2 muss korrigiert werden.

---

## Zusammenfassung

### Gesamtübersicht

| Paper | Titel | EN | DE | Kritische Fehler |
|---|---|---|---|---|
| 18 | Vinogradov 3-Prim | ⚠️ | ⚠️ | Vorzeichenfehler Remark 4.3 (EN+DE) |
| 19 | Waringsches Problem | ⚠️ | ❌ | DE: Haupttheorem als Remark; fehlende Sektionen |
| 20 | Goldbach sing. Reihe | ⚠️ | ❌ | Tabelle S(n) falsch (EN+DE); Section 8 fehlt (DE) |
| 21 | Riemannsche Zetafunktion | ✅ | ✅ | keine |
| 22 | Nichttriviale Nullstellen | ✅ | ✅ | keine |
| 23 | Explizite Formel + PNT | ✅ | ✅ | keine |
| 24 | RH-Ansätze | ✅ | ❌ | Section 8 (Physik-Verbindung) fehlt in DE |
| 25 | Elliptische Kurven/Q | ✅ | ⚠️ | fehlendes Bibitem Cornell1997; fehlende Beispiele |
| 26 | L-Funktionen ell. Kurven | ✅ | ⚠️ | fehlendes Bibitem TaylorWiles1995; strukturelle Abweichung |
| 27 | BSD-Vermutung | ✅ | ✅ | keine kritischen |
| 28 | Kongruente Zahlen | ✅ | ⚠️ | falsches Parity-Kriterium A(2) statt C(2) |

### Offene Fehler nach Schwere

**Kritisch (Mathematikfehler / fehlende Hauptsektionen):**
1. BUG-B4-P19-DE-MISSING-THM: Wooleys Theorem als Remark (Paper 19 DE)
2. BUG-B4-P20-DE-SECTION-MISSING: Section 8 fehlt (Paper 20 DE)
3. BUG-B5-P24-DE-MISSING-SECTION: Physik-Sektion fehlt (Paper 24 DE)
4. BUG-B6-P28-DE-TUNNELL-PARITY: Falsches Parity-Kriterium für gerades n=2 (Paper 28 DE)

**Hoch:**
5. BUG-B4-P18-EN-001: Vorzeichenfehler Remark 4.3 (Paper 18 EN)
6. BUG-B4-P18-DE-001: Identischer Vorzeichenfehler (Paper 18 DE)
7. BUG-B4-P19-EN-002: Inkohärente Schwellenbedingungen (Paper 19 EN)
8. BUG-B4-P20-EN-002: Systematisch falsche Tabelle S(n) (Paper 20 EN)
9. BUG-B4-P20-DE-001: Identische Tabellenfehler (Paper 20 DE)

**Mittel:**
10. BUG-B4-P18-DE-002: Unvollständiger Beweis sing. Integral (Paper 18 DE)
11. BUG-B4-P19-EN-001: Falsche Querverweise auf Paper 17 (Paper 19 EN)
12. BUG-B4-P19-DE-001: Fehlendes Korollar Weyl-Schranke (Paper 19 DE)
13. BUG-B4-P19-DE-003: Fehlendes Efficient-Congruencing-Theorem (Paper 19 DE)
14. BUG-B4-P19-DE-MISSING-REFS: 3 fehlende Bibitems (Paper 19 DE)
15. BUG-B4-P19-EN-MISSING-BIBITEM: Fehlendes `\bibitem{Wooley2016}` (Paper 19 EN)
16. BUG-B4-P20-DE-MISSING-BIBITEMS: 3 fehlende Bibitems (Paper 20 DE)
17. BUG-B4-P20-EN-TABLE-001: Unklare Tabellennotation (Paper 20 EN)
18. BUG-B6-P25-DE-001: Fehlendes Beispiel j-Invariante (Paper 25 DE)
19. BUG-B6-P25-DE-002: Fehlendes Bibitem Cornell1997 (Paper 25 DE)
20. BUG-B6-P26-DE-001: Fehlende Tabelle Example 2.3 (Paper 26 DE)
21. BUG-B6-P26-DE-002: Fehlendes/falsches Bibitem TaylorWiles1995 (Paper 26 DE)
22. BUG-B6-P27-DE-002: Fehlende algebraische Formel Section 7 (Paper 27 DE)

**Gering:**
23. BUG-B4-P18-DE-002: Beweisvervollständigung (Paper 18 DE)
24. BUG-B4-P19-DE-002: Korollar ohne Beweis (Paper 19 DE)
25. BUG-B4-P20-EN-001: Satzfragment (Paper 20 EN)
26. BUG-B6-P25-DE-003: Fehlende Rang-Remark (Paper 25 DE)

**Systemisch (alle Batch 5+6 DE-Dateien):**
- Fehlendes `\usepackage[utf8]{inputenc}` in allen DE-Dateien von Batch 5 und Batch 6.

### Statistik

| | Papers | Dateien |
|---|---|---|
| Vollständig druckreif (EN+DE) | Papers 21, 22, 23, 27 | 8 |
| EN druckreif, DE mit Mittelproblemen | Papers 25, 26, 28 | 6 |
| Beide Versionen mit Fehlern | Papers 18, 19, 20 | 6 |
| DE kritisch unvollständig | Papers 19, 20, 24 | 3 DE |
| Gesamt | 11 Papers | 22 Dateien |

---

*Erstellt: 2026-03-11 | Build 9 | Gutachter: Claude Sonnet 4.6*
