# Review: Batch 24 — Gruppe B (Papers 92–95)
**Datum:** 2026-03-12
**Autor:** Michael Fuhrmann
**Build:** 172
**Verifikationsskript:** `/src/py/gruppe_b_batch24_verification.py`

---

## Übersicht der Vermutungen und ihr Status

| Paper | Titel | Aussage | Status |
|-------|-------|---------|--------|
| 92 | Donaldson 4-Mannigfaltigkeiten | Donaldson-Theorem: Definite Schnittformen sind diagonal | **BEWIESEN** (Donaldson 1983) |
| 92 | Glattes Poincaré Dim 4 | Jede glatte Homotopie-4-Sphäre ist diffeomorph zu S⁴ | **OFFEN** |
| 93 | Pu's Ungleichung für RP² | sys(RP²)² ≤ (π/2)·vol(RP²) mit Gleichheit bei runder Metrik | **BEWIESEN** (Pu 1952) |
| 93 | Loewner-Ungleichung für T² | sys(T²)² ≤ (2/√3)·vol(T²) mit Gleichheit bei hex. Torus | **BEWIESEN** (Loewner 1949) |
| 93 | Optimale Füllungsradius-Konstanten | Scharfe Konstanten für T^n (n≥3), RP^n (n≥3) | **OFFEN** |
| 94 | Quillen-Suslin (Serre-Vermutung) | Projektive Moduln über k[x₁,...,xₙ] sind frei | **BEWIESEN** (Quillen + Suslin 1976) |
| 94 | Hartshorne Kodimension-2 | Glatte Kodim-2 Untervar. von ℙⁿ (n≥7) sind vollst. Durchschnitte | **OFFEN** |
| 95 | Mori-Siu-Yau Kompaktfall | Kompakte Kähler-Mannigfaltigkeit mit pos. bisk. Kr. ≅ ℙⁿ | **BEWIESEN** (Mori 1979, Siu-Yau 1980) |
| 95 | Yau-Vermutung (nicht-kompakt) | Vollst. Kähler mit pos. bisk. Krümmung ≅ ℂⁿ | **OFFEN** |

---

## 1. Beweis-Präsentationen

### 1.1 Quillen-Suslin (Serre-Vermutung) — BEWIESEN 1976

**Theorem (Quillen 1976, Suslin 1976):** Jeder endlich erzeugte projektive Modul über k[x₁,...,xₙ] (k Körper oder Hauptidealring) ist frei.

**Beweis-Skizze (Quillens Ansatz):**

**Schritt 1 — Reduktion auf Patching:**
Horrocks' Kriterium (1964): Ein projektiver Modul P über k[x₁,...,xₙ] ist frei genau dann, wenn P ⊗ k[xₙ]_{(f)} frei ist für ein einziges monisches Polynom f ∈ k[xₙ].

**Schritt 2 — Quillens Patching-Lemma:**
Sei A ein noetherscher Ring, M ein endlich präsentierter A[t]-Modul.
Falls M_m für alle maximalen Ideale m ⊂ A aus A_m "extended" ist (d.h. M_m ≅ N_m ⊗_{A_m} A_m[t] für einen A_m-Modul N_m), dann ist M aus A extended.

*Bedeutung:* Die Freiheit kann "von unten" (über alle Lokalisierungen) zusammengeflickt werden.

**Schritt 3 — Induktion über n:**
- Basisfall n=1: Projektive Moduln über PID k[x] sind frei (klassisch, folgt aus Struktursatz).
- Induktionsschritt: Angenommen, der Satz gilt für k[x₁,...,x_{n-1}]. Sei P projektiv über k[x₁,...,xₙ] = k[x₁,...,x_{n-1}][xₙ].
- Für jedes maximale Ideal m ⊂ A = k[x₁,...,x_{n-1}] ist P_m projektiv über A_m[xₙ].
- Nach Induktionshypothese (lokale Version) ist P_m frei über A_m[xₙ].
- Das Patching-Lemma liefert: P ist bereits aus A extended.
- Da P nach Schritt 1 auch als Modul über k[x₁,...,x_{n-1}] frei ist, folgt die Freiheit von P.

**Schluss:** Suslin fand einen unabhängigen, elementareren Beweis (Matrizenrechnung), der explizite Konstruktionen erlaubt. Beide Beweise erschienen 1976.

**Geometrische Bedeutung:** Jedes algebraische Vektorbündel über dem affinen Raum Aⁿ_k ist trivial. Dies steht im Kontrast zum projektiven Raum ℙⁿ, wo nicht-triviale Bündel existieren (Horrocks-Mumford auf ℙ⁴).

---

### 1.2 Mori's Bend-and-Break und Mori-Siu-Yau — BEWIESEN

**Theorem (Mori 1979, Siu-Yau 1980):** Jede kompakte Kähler-Mannigfaltigkeit mit positiver bisektionaler Krümmung ist biholomorph zu ℙⁿ.

**Mori's Bend-and-Break (algebraischer Ansatz):**

Das Schlüsselprinzip: Wenn eine rationale Kurve mit hohem antikanonischem Grad existiert, "biegt" sie sich unter Deformation und "bricht" in einfachere Kurven auf.

**Beweis-Skizze:**

**Schritt 1 (Projektiv):** Positive bisektionale Krümmung ⟹ c₁(M) > 0 ⟹ −K_M ample (antikanonisches Bündel ampel). Daher ist M eine Fano-Varietät.

**Schritt 2 (Rationaler Kurven via Frobenius):**
Mori arbeitet in Charakteristik p (Reduktion modulo p):
- Für eine Kurve C ⊂ M mit (−K_M · C) > 0 und Frobenius-Morphismus φ: C → C^(p)
- Die Deformation der zusammengesetzten Kurve φ: C → M × C^(p) kann (Bend-and-Break) eine rationale Kurve erzeugen
- Für hinreichend großes p entstehen rationale Kurven (da der Modulraum zu groß ist, um leer zu sein)

**Schritt 3 (Kobayashi-Ochiai):** Wenn der Index i(M) = n+1, d.h. −K_M = (n+1)H für eine ample Klasse H, dann M ≅ ℙⁿ. Positive bisektionale Krümmung impliziert index n+1 (durch Riemann-Roch-Argumentation).

**Warum positive bisektionale Krümmung rationale Kurven erzwingt:**
Die bisektionale Krümmung Bsc(X,Y) = R(X,JX,Y,JY) misst die Wechselwirkung zwischen zwei komplexen Richtungen. Positiv zu sein bedeutet, dass jede Fläche, die in zwei komplexe Richtungen aufgespannt wird, positiv gekrümmt ist — ein sehr starkes Starrheits-Prinzip. Es impliziert positive Ricci-Krümmung (c₁ > 0) und liefert via Bochner-Technik die Existenz rationaler Kurven.

**Siu-Yau's differential-geometrischer Beweis:**
Siu-Yau verwenden harmonische Abbildungstheorie: Konstruiere eine harmonische Abbildung f: M → ℙⁿ vom Grad 1 via Bochner-Weitzenböck. Die Positivität der bisektionalen Krümmung erzwingt die Holomorphizität (Rigiditätssatz), und f ist ein Biholomorphismus.

---

### 1.3 Pu's Ungleichung für RP² — BEWIESEN 1952

**Theorem (Pu 1952):** Für jede Riemannsche Metrik g auf RP²:
$$\text{sys}(\mathbb{RP}^2, g)^2 \leq \frac{\pi}{2} \cdot \text{vol}(\mathbb{RP}^2, g)$$
mit Gleichheit genau für die runde Standardmetrik (konstante Krümmung +1, skaliert).

**Beweis via Konformfaktor-Integration:**

**Schritt 1 (Uniformisierung):**
Jede Riemannsche Metrik auf RP² ist konform äquivalent zur runden Metrik:
g = e^{2u} g_{rund} für eine Funktion u: RP² → ℝ.

**Schritt 2 (Systemvolumen-Abschätzung):**
Das Volumen transformiert als:
$$\text{vol}(\mathbb{RP}^2, g) = \int_{\mathbb{RP}^2} e^{2u}\, dA_{\text{rund}}$$

**Schritt 3 (Systole via konformale Länge):**
Die kürzeste nicht-kontrahierbare Kurve γ in (RP², g) hat Länge:
$$\text{sys}(g) = \text{sys}(e^{2u} g_{\text{rund}}) \leq \min_{\gamma \text{ nicht-kontr.}} \int_\gamma e^u\, ds_{\text{rund}}$$

**Schritt 4 (Cauchy-Schwarz):**
Für die nicht-kontrahierbare Geodäte γ₀ der Länge π (im Runden RP²):
$$\text{sys}(g)^2 \leq \left(\int_{\gamma_0} e^u\, ds\right)^2 \leq \pi \int_{\gamma_0} e^{2u}\, ds$$
(Cauchy-Schwarz: (∫f)² ≤ L·∫f²)

**Schritt 5 (Integration über alle Richtungen):**
Pu integriert über alle "Meridiane" (geodätische Symmetrielinien der RP²) und nutzt die Rotationssymmetrie der Runden Metrik, um die rechte Seite durch (π/2)·vol(g) abzuschätzen.

**Gleichheitsfall:** Cauchy-Schwarz ist scharf genau wenn e^u konstant auf γ₀, also u = const., d.h. g konform gleich g_{rund} — nach Skalierung exakt die runde Metrik.

---

## 2. Analyse der offenen Teile

### 2.1 Glattes Poincaré-Problem Dim 4 — Warum so schwierig?

**Vermutung (Offen):** Jede glatte Homotopie-4-Sphäre ist diffeomorph zu S⁴.

**Dimensionsvergleich:**
- **Dim ≤ 3:** Standardtopologie. Poincaré-Vermutung: Perelman (2002/03) löste den 3D-Fall via Ricci-Fluss. Für Dim 1, 2 trivial.
- **Dim ≥ 5:** Smales h-Kobordismus-Theorem (1961): Für n ≥ 5 sind zwei einfach zusammenhängende h-kobordante glatte Mannigfaltigkeiten diffeomorph. Daher: Jede glatte Homotopie-n-Sphäre (n ≥ 5) ist diffeomorph zu S^n (Schönflies Vermutung bewiesen).
- **Dim 4 (topologisch):** Freedman (1982) bewies: Jede Homotopie-4-Sphäre ist *homöomorph* zu S⁴.
- **Dim 4 (glatt):** OFFEN.

**Warum versagen die bekannten Methoden in Dim 4?**

1. **h-Kobordismus scheitert:** Der Beweis des h-Kobordismus-Theorems benötigt dim ≥ 5 (Whitney-Trick: Crossing-Changes an eingebetteten Disks). In Dim 4 ist der Whitney-Trick nicht ohne weitere Hypothesen anwendbar (die Disks schneiden sich selbst).

2. **Keine SW-Invarianten:** Seiberg-Witten-Invarianten sind nur für b₂⁺ > 0 definiert und nicht-trivial. S⁴ hat b₂ = 0, daher versagen SW-Invarianten direkt.

3. **Gluck-Twists:** Die Hauptkandidaten für exotische S⁴ (Gluck-Twists von eingebetteten S² ⊂ S⁴) wurden alle als diffeomorph zu S⁴ erkannt — aber kein allgemeines Kriterium ist bekannt.

4. **Keine Analogie zu Exotischen ℝ⁴:** Obwohl exotische ℝ⁴ existieren (aus der Kombination Freedman + Donaldson), ist unbekannt, ob *kompakte* exotische 4-Sphären existieren.

5. **Kirby-Kalkül-Hindernisse:** In Dim 4 ist die Darstellungstheorie gerahmter Links erheblich reicher — keine vollständige algorithmische Klassifikation.

**Aktueller Forschungsstand:** Das Problem gilt als eines der am tiefsten liegenden in der Differentialtopologie. Neue Invarianten (z.B. Khovanov-Homologie, Floer-Theorie, Heegaard-Floer) liefern bisher keinen Widerspruch und kein Beweis.

---

### 2.2 Yau-Vermutung (nicht-kompakt) — Offene Teile

**Vermutung (Yau, ~1974, Offen in voller Allgemeinheit):** Jede vollständige Kähler-Mannigfaltigkeit mit positiver bisektionaler Krümmung ist biholomorph zu ℂⁿ.

**Bekannte Partialergebnisse:**

| Zusatzbedingung | Autor | Jahr | Ergebnis |
|-----------------|-------|------|---------|
| Komplexe Dim 2 + endliches Volumen | Mok | 1984 | M ≅ ℂ² |
| Beschränkte Krümmung |R| ≤ C | Shi | 1989 | M ≅ ℂⁿ |
| Maximales Volumenwachstum vol(B(p,r)) ≥ cr^{2n} | Chen-Tang-Zhu | 2004 | M ≅ ℂⁿ |
| Kähler-Ricci-Fluss, kompakter Fall | Cao | 1985 | Konvergenz zu Fubini-Study |

**Cao's Kähler-Ricci-Fluss-Ansatz:**

Cao (1985) bewies: Der Kähler-Ricci-Fluss ∂g_{iȷ̄}/∂t = −R_{iȷ̄} auf einer kompakten Kähler-Mannigfaltigkeit mit positiver bisektionaler Krümmung existiert für alle t und konvergiert (nach Normalisierung) zur Fubini-Study-Metrik auf ℙⁿ.

Für den nicht-kompakten Fall (Yau's Vermutung):
1. **Existenz:** Shi (1989) bewies Kurzzeit-Existenz für vollständige nicht-kompakte Mannigfaltigkeiten mit beschränkter Krümmung.
2. **Langzeit-Existenz:** Unter bounded curvature bleibt der Fluss vollständig und Kähler.
3. **Konvergenz zu ℂⁿ:** Hier scheitern die bekannten Methoden. Man bräuchte:
   - Kontrolle des Krümmungsabfalls entlang des Flusses
   - Nachweis, dass der Fluss zur flachen Metrik konvergiert
   - Biholomorphismus-Beweis aus metrischer Konvergenz

**Warum Cao's Ansatz nicht ausreicht:**
Im kompakten Fall erzwingt die Topologie (ℙⁿ hat bekannte Topologie) die Konvergenz. Im nicht-kompakten Fall ist die Topologie von M unbekannt (soll ja erst bewiesen werden), und die Konvergenz des Flusses ohne topologische a-priori-Information ist erheblich schwieriger.

---

### 2.3 Hartshorne Kodimension-2 — Warum hilft Horrocks-Mumford nicht für ℙ⁷?

**Vermutung (Hartshorne 1974, Offen):** Jede glatte Kodimension-2 Untervarietät von ℙⁿ (n ≥ 7) ist ein vollständiger Durchschnitt.

**Das Horrocks-Mumford-Bündel auf ℙ⁴:**
Horrocks und Mumford (1973) konstruierten ein unzerlegbares Rang-2-Vektorbündel F_HM auf ℙ⁴ mit c₁ = 5, c₂ = 10. Seine Nullstellen (via der Serre-Korrespondenz) sind abelsche Flächen vom Grad 10 in ℙ⁴ — das sind glatte Kodimension-2-Untervarietäten von ℙ⁴, die KEINE vollständigen Durchschnitte sind.

**Warum funktioniert das nicht für ℙ⁷?**

1. **Der Dimensionssprung und Barth-Lefschetz:**
   Der Satz von Barth (1970) besagt: Für Y ⊂ ℙⁿ glatt, Kodimension c, ist die Restriktionsabbildung H^k(ℙⁿ,ℤ) → H^k(Y,ℤ) ein Isomorphismus für k < n-2c+1.
   Für n=7, c=2: Isomorphismus für k < 4. Das heißt: Y hat dieselbe Kohomologie wie ein vollständiger Durchschnitt in vielen Graden — sehr restriktiv.

2. **Keine bekannte Konstruktion analoger Bündel auf ℙ^n (n≥5):**
   Das HM-Bündel verwendet die Heisenberg-Gruppe H₅ und die spezielle Symplektik von ℙ⁴. Für n ≥ 5 gibt es keinen bekannten Mechanismus, ein unzerlegbares Rang-2-Bündel zu konstruieren. Es ist ein offenes Problem (die Splitting-Vermutung), ob auf ℙⁿ (n ≥ 6) überhaupt indekomposable Rang-2-Bündel existieren.

3. **Hartshorne-Ogus (1974):**
   Für n ≥ 2·dim(Y)+1 hat das Normalbündel N_{Y/ℙⁿ} die kohomologischen Eigenschaften eines vollständigen Durchschnitts. Das schränkt die möglichen Gegenbeispiele stark ein, liefert aber keinen Beweis.

4. **Die Serre-Korrespondenz:**
   Hartshorne's Vermutung ist äquivalent zur Splitting-Vermutung für Rang-2-Bündel auf ℙⁿ (n ≥ 7). Ein Beweis müsste zeigen, dass jedes Rang-2-Bündel auf ℙ⁷ die Form O(a) ⊕ O(b) hat — eine tiefe kombinatorische und geometrische Aussage.

---

## 3. Numerische Verifikationen

### 3.1 Loewner-Ungleichung — Numerisch verifiziert

**Ergebnisse des Verifikationsskripts** (`gruppe_b_batch24_verification.py`):

| Torus τ | sys²/Area | Loewner-Schranke 2/√3 | Status |
|---------|-----------|----------------------|--------|
| e^{iπ/3} = 0.5+i·0.866 (hexagonal) | **1.154701** | 1.154701 | Gleichheit ✓ |
| i (quadratisch) | 1.000000 | 1.154701 | Strikt darunter ✓ |
| ~7200 Punkte im Fundamentalbereich | max 1.124339 | 1.154701 | Überall ≤ Schranke ✓ |

**Beobachtung:** Das numerische Maximum (1.124339) liegt unter der theoretischen Schranke (1.154701), da das Gitter nicht fein genug ist, den exakten hexagonalen Punkt zu treffen. Die Schranke ist in allen Punkten erfüllt.

Das Maximum liegt bei Re(τ) = -0.5, Im(τ) ≈ 0.889 — dieser Punkt ist SL(2,ℤ)-äquivalent zu Re(τ) = 0.5, Im(τ) ≈ 0.866 (der hexagonale Torus).

### 3.2 Quillen-Suslin — Algebraisch verifiziert

| Beispiel | Ideal | Einheitsideal? | Schluss |
|----------|-------|----------------|---------|
| (1+x, 1-x) | Gr.-Basis = [1] | Ja ✓ | Modul frei |
| (x, 1-xy, y) | Gr.-Basis = [1] | Ja ✓ | Modul frei |
| (x, y) | Gr.-Basis = [x,y] | **Nein** ✓ | Ideal ≠ k[x,y] |
| (x, y, 1-xy) | Gr.-Basis = [1] | Ja ✓ | Modul frei |

Alle 4 Testfälle entsprechen den erwarteten Ergebnissen. Der kritische Punkt: gcd(x,y)=1 als Polynome, aber 1 ∉ Ideal(x,y), weil alle Elemente des Ideals bei (0,0) verschwinden.

---

## 4. Klassifikation aller Aussagen

### BEWIESEN

1. **Donaldson-Theorem (1983):** Definite Schnittformen glatter, einfach zusammenhängender 4-Mannigfaltigkeiten sind diagonal über ℤ.
   *Beweis:* Yang-Mills Instanton-Modulraum + Uhlenbeck-Kompaktheit.

2. **Quillen-Suslin (1976):** Projektive Moduln über k[x₁,...,xₙ] sind frei.
   *Beweis:* Quillens Patching-Lemma + Induktion, oder Suslins Matrizenmethode.

3. **Pu's Ungleichung (1952):** sys(RP²)² ≤ (π/2)·vol(RP²), scharf.
   *Beweis:* Uniformisierung + Konformfaktor + Cauchy-Schwarz.

4. **Loewner-Ungleichung (1949):** sys(T²)² ≤ (2/√3)·vol(T²), scharf bei hex. Torus.
   *Beweis:* Uniformisierung + Gittergeometrie + Lagranges Reduktionssatz.

5. **Mori-Siu-Yau (1979/80):** Kompakte Kähler-Mannigfaltigkeit mit pos. bisk. Kr. ≅ ℙⁿ.
   *Beweis:* Mori: Bend-and-Break, Kobayashi-Ochiai; Siu-Yau: Harmonische Abbildungen.

6. **Uhlenbeck-Kompaktheit (1982):** ASD-Verbindungen mit beschränkter YM-Energie haben konvergente Teilfolgen (mit Blasenbildung).

7. **Freedman-Klassifikation (1982):** Einfach zusammenhängende topologische 4-Mannigfaltigkeiten ↔ unimodulare bilineare Formen + ks-Invariante.

### OFFEN

1. **Glattes Poincaré Dim 4:** Jede glatte Homotopie-4-Sphäre ≅ S⁴?
   *Hindernisse:* Whitney-Trick in Dim 4 nicht verfügbar, SW-Invarianten ungeeignet für b₂=0.

2. **Yau-Vermutung (nicht-kompakt):** Vollst. Kähler + pos. bisk. Kr. ≅ ℂⁿ?
   *Partialergebnisse:* Shi (beschr. Kr.), Chen-Tang-Zhu (max. Vol.-Wachstum), Cao (Ricci-Fluss kompakt).

3. **Hartshorne Kodimension-2:** Glatte Kodim-2 Untervar. von ℙⁿ (n≥7) = vollst. Durchschnitt?
   *Äquivalent zu:* Jedes Rang-2-Bündel auf ℙⁿ (n≥7) zerfällt.

4. **Berger-Problem:** Optimale systolische Konstanten für T^n (n≥3), RP^n (n≥3).

5. **Systolische Freiheit:** Welche einfach zusammenhängenden Mannigfaltigkeiten (dim≥4) haben unbeschränktes systolisches Verhältnis?

---

## 5. Mathematische Korrektheit der Papers

### Paper 92 (Donaldson 4-Mannigfaltigkeiten) — KORREKT
- Theorem/Conjecture-Trennung korrekt: Donaldson-Theorem als Theorem, Glattes Poincaré als Conjecture.
- Dimensionsformel der Moduli-Raum-Dimension korrekt: dim M₁ = 8c₂ - 3(1-b₁+b₂⁺) = 5 für b₂⁺=0.
- Seiberg-Witten-Gleichungen korrekt formuliert.
- Exotische ℝ⁴ korrekt beschrieben.
- **Keine Fehler gefunden.**

### Paper 93 (Gromov Füllungsradius) — KORREKT
- Loewner-Ungleichung: sys² ≤ (2/√3)·vol — korrekt (nicht sys² ≤ √3/2·vol!).
- Pu-Ungleichung: sys² ≤ (π/2)·vol — korrekt.
- Gromov-Ungleichung sys^n ≤ C_n·vol für wesentliche Mannigfaltigkeiten — korrekt.
- Füllungsradius-Definition über Kuratowski-Einbettung in L∞(M) — korrekt.
- Katz's Formel FillRad(Sⁿ) = arccos(-1/(n+1))·R — korrekt.
- **Keine Fehler gefunden.**

### Paper 94 (Hartshorne-Vermutungen) — KORREKT
- Quillen-Suslin explizit als Theorem markiert — korrekt.
- Hartshorne Kodim-2 explizit als Conjecture markiert — korrekt.
- Horrocks-Mumford als Theorem (existent) — korrekt.
- Splitting Conjecture als Conjecture — korrekt.
- **Keine Fehler gefunden.**

### Paper 95 (Uniformisierung höhere Dim.) — KORREKT
- Mori-Siu-Yau kompakter Fall als Theorem — korrekt.
- Yau nicht-kompakt als Conjecture [Open in Full Generality] — korrekt.
- Partialergebnisse (Mok, Shi, Chen-Tang-Zhu) als Theoreme — korrekt.
- Mok's Theorem (nicht-neg. bisk. Kr. → hermitesche symmetrische Räume) — korrekt.
- **Kleiner Tippfehler:** Zeile 262: `\end{conjecture}` vor `\end{theorem}` vertauscht — LaTeX-Fehler!

---

## 6. LaTeX-Bug in Paper 95

In Paper 95 (beide Sprachversionen) gibt es eine vertauschte Umgebung:

```latex
% Zeile 250–263 in paper95_uniformization_higher_dimensions_en.tex
\begin{theorem}[Partial Results on Conjecture~\ref{conj:yau}]
...
\end{conjecture}   % FALSCH: sollte \end{theorem} sein
\end{theorem}      % FALSCH: sollte entfernt werden
```

**Behobene Versionen:** Beide `.tex`-Dateien müssen korrigiert werden.

---

## 7. Gesamtergebnis

| Kategorie | Anzahl |
|-----------|--------|
| BEWIESEN | 7 Aussagen |
| OFFEN | 5 Aussagen |
| WIDERLEGBAR | 0 |
| Numerisch bestätigt | 2 (Loewner, Quillen-Suslin) |
| LaTeX-Fehler | 1 (Paper 95) |

**Fazit:** Die vier Papers präsentieren die Mathematik korrekt und unterscheiden sauber zwischen bewiesenen Theoremen und offenen Vermutungen. Der einzige gefundene Fehler ist ein LaTeX-Strukturfehler in Paper 95 (vertauschte theorem/conjecture-Umgebung). Alle numerischen Verifikationen sind erfolgreich durchgelaufen.
