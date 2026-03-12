# Mathematisches Audit: Batch 5 (Riemann Zeta) & Batch 6 (Elliptische Kurven / BSD)
## Gutachter: Claude Sonnet 4.6 | Datum: 2026-03-12 | Build: 10

---

## BATCH 5 — RIEMANN ZETA-FUNKTION (Papers 21–24)

---

### Paper 21 (EN): The Riemann Zeta Function — Definition, Analytic Continuation, and the Functional Equation

**Hauptbehauptung:** Entwicklung der Zeta-Funktion von Grund auf: Dirichlet-Reihe, Euler-Produkt, analytische Fortsetzung auf C\{1}, Funktionalgleichung, triviale Nullstellen, Riemann–Siegel-Formel, Euler-Werte ζ(2k).

**Kritische Frage:** RH wird durchgehend als OFFENES PROBLEM korrekt behandelt. Explizit: *"the Riemann Hypothesis (still open)"* (Abstract), *"remains the most celebrated open problem in mathematics"* (Schluss).

**Mathematische Korrektheit:** ✅

**Prüfung der Kernformeln:**

1. **Funktionalgleichung** (eq. functional):
   `ζ(s) = 2^s · π^(s-1) · sin(πs/2) · Γ(1-s) · ζ(1-s)` ✅
   Äquivalente ξ-Form: `ξ(s) = ½ s(s-1) π^(-s/2) Γ(s/2) ζ(s)`, `ξ(s) = ξ(1-s)` ✅

2. **Residuum bei s=1:** Laurent-Expansion `ζ(s) = 1/(s-1) + γ + O(s-1)` ✅

3. **Triviale Nullstellen:** Beweis via sin(-mπ) = 0 korrekt. Einfachheit der Nullstellen korrekt begründet. ✅

4. **Gamma-Funktion:** Reflexionsformel `Γ(s)Γ(1-s) = π/sin(πs)` ✅. Pole bei s = 0, -1, -2, ... ✅

5. **Riemann–Siegel-Formel:**
   - θ(t) = arg[Γ(1/4 + it/2)] - (t/2)·ln π ✅
   - Asymptotik: θ(t) = (t/2)ln(t/2π) - t/2 - π/8 + 1/(48t) - 7/(5760t³) + O(t^-5) ✅
   - Z(t) = e^(iθ(t)) · ζ(1/2 + it) ✅
   - N(t) = ⌊√(t/2π)⌋ ✅

6. **Euler-Werte:**
   `ζ(2k) = (-1)^(k+1)(2π)^(2k) B_{2k} / (2(2k)!)` ✅
   Spezialwerte: ζ(2) = π²/6, ζ(4) = π⁴/90, ζ(6) = π⁶/945 ✅

7. **ζ(-1) = -1/12, ζ(-3) = 1/120, ζ(0) = -1/2** ✅

8. **Abschätzung bei s=0:** Behauptung ζ(0) = -1/2 ≠ 0 korrekt, Beweis dass s=0 keine Nullstelle ist, korrekt. ✅

**Fehler:** Keine mathematischen Fehler gefunden.

**Urteil: DRUCKREIF**

---

### Paper 21 (DE): Die Riemann'sche Zeta-Funktion

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Vollständige und korrekte Übersetzung. Alle Formeln identisch. Keine inhaltlichen Abweichungen. ✅

**Einzige formale Anmerkung:** Im DE-Paper fehlt `\usepackage[utf8]{inputenc}` — NEIN, tatsächlich vorhanden in Zeile 5. ✅

**Fehler:** Keine.

**Urteil: DRUCKREIF**

---

### Paper 22 (EN): The Non-Trivial Zeros of the Riemann Zeta Function — Counting, Distribution, and the GUE Conjecture

**Hauptbehauptung:** Nullstellen-Zählformel N(T), erste nichttriviale Nullstelle, Gram-Punkte, Montgomery–Odlyzko GUE-Statistik, numerische Verifikation bis 10^13, GRH.

**Kritische Frage:** RH wird durchgehend als OFFENES PROBLEM behandelt. *"the Riemann Hypothesis remains one of the greatest unsolved problems in mathematics."* ✅

**Mathematische Korrektheit:** ✅

**Prüfung der Kernformeln:**

1. **Riemann–von Mangoldt-Formel:**
   `N(T) = (T/2π)log(T/2π) - T/2π + 7/8 + O(log T)` ✅
   (Im Abstract ist O(log T) ohne die Konstante 7/8 — das ist eine akzeptable Kurzform des Hauptterms, da 7/8 zum O(log T)-Term gehört; die genaue Formel mit 7/8 steht korrekt in Theorem 2.1.)

2. **Erste Nullstelle:** γ₁ = 14.134725141734693... ✅ (anerkannter Wert, LMFDB-bestätigt)

3. **Zehn erste Nullstellen:** Alle Werte mit korrekten Stellen. ✅

4. **Normalisierter Nullstellenabstand:** δₙ = (γₙ₊₁ - γₙ)·log(γₙ)/(2π) ✅

5. **GUE-Paarkorrelation:** r₂(s) = 1 - (sin(πs)/πs)² ✅

6. **Montgomery-Formel:** Korrekter bedingter Satz (RH-Voraussetzung explizit), Einschränkung auf Testfunktionen mit Träger in (-1,1) in Remark korrekt erwähnt. ✅

7. **Gram-Punkte:** θ(gₙ) = nπ ✅. Erste Ausnahme bei n=126 ✅.

8. **GRH für L(s,χ):** Schranke p ≤ C·q²(log q)² korrekt als GRH-Konsequenz. ✅

**Fehler:** Keine mathematischen Fehler.

**Urteil: DRUCKREIF**

---

### Paper 22 (DE): Die nicht-trivialen Nullstellen der Riemann'schen Zeta-Funktion

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Vollständige inhaltliche Übereinstimmung. Alle Formeln identisch. ✅

**Fehler:** Keine.

**Urteil: DRUCKREIF**

---

### Paper 23 (EN): Riemann's Explicit Formula and the Prime Number Theorem

**Hauptbehauptung:** Explizite Formel für ψ(x) und π(x), Primzahlsatz via nullstellenfreier Bereich Re(s)=1, Fehlerterm unter RH, Schoenfeld-Schranke, GRH-Anwendung.

**Kritische Frage:** RH wird nur als Annahme verwendet, nie als bewiesen behandelt. ✅ PNT wird korrekt als *bewiesen* dargestellt (Hadamard/de la Vallée-Poussin 1896). ✅

**Mathematische Korrektheit:** ✅

**Prüfung der Kernformeln:**

1. **Explizite Formel für ψ(x):**
   `ψ(x) = x - Σ_ρ x^ρ/ρ - ζ'(0)/ζ(0) - ½ log(1-x^{-2})`
   mit `ζ'(0)/ζ(0) = log(2π)` ✅

   **WICHTIGER HINWEIS:** Die Standardformel bei Davenport/Titchmarsh lautet
   `ψ(x) = x - Σ_ρ x^ρ/ρ - log(2π) - ½ log(1-x^{-2})`
   Das Paper schreibt `-ζ'(0)/ζ(0)`, und gibt an `ζ'(0)/ζ(0) = log(2π)`.
   Daher lautet der Term: `-log(2π)`, was korrekt ist. ✅

   Der Beweis-Kommentar erklärt den s=0-Beitrag sehr sorgfältig (BUG-B5-P23-EN/DE-001 Kommentar im Code): ζ(0) = -1/2 ≠ 0, daher hat -ζ'/ζ keinen Pol bei s=0; der 1/s-Faktor liefert das Residuum. Das ist mathematisch korrekt und präzise. ✅

2. **Riemann-Explizitformel für π(x):**
   `π(x) = Li(x) - Σ_ρ Li(x^ρ) - log 2 + ∫_x^∞ dt/(t(t²-1)log t)` ✅
   mit Li(x) = ∫₂^x dt/log t ✅

3. **Mertens-Argument (keine Nullstellen auf Re(s)=1):**
   `log|ζ³(σ)·ζ⁴(σ+it)·ζ(σ+2it)| ≥ 0` via `3+4cosθ+cos2θ = 2(1+cosθ)² ≥ 0` ✅
   Der Widerspruchsbeweis ist korrekt formuliert. ✅

4. **Fehlerterm:** `ψ(x) = x + O(x^(β₀)(log x)²)` unter β₀-Schranke ✅
   Unter RH: `ψ(x) = x + O(√x log²x)` ✅

5. **Schoenfeld-Schranke:** `|π(x)-Li(x)| < (1/8π)√x log x` für x ≥ 2657 ✅

6. **Logarithmische Ableitung:** `-ζ'(s)/ζ(s) = Σ Λ(n)/n^s` für Re(s)>1 ✅

7. **Perron-Formel:** Korrekte Darstellung als Bromwich-Integral. ✅

8. **GRH-Fehlerterm:** `π(x;q,a) = Li(x)/φ(q) + O(√x log(qx))` unter GRH ✅

**Fehler:** Keine mathematischen Fehler.

**Urteil: DRUCKREIF**

---

### Paper 23 (DE): Riemanns explizite Formel und der Primzahlsatz

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Vollständige und korrekte Übersetzung. Alle Formeln identisch. ✅
Zusätzlich enthält die DE-Version einen ausführlicheren Abstract, der GRH-Anwendungen expliziter beschreibt. Das ist eine sinnvolle Erweiterung. ✅

**Fehler:** Keine.

**Urteil: DRUCKREIF**

---

### Paper 24 (EN): Approaches to the Riemann Hypothesis — Hilbert–Pólya, Berry–Keating, Random Matrices, and Open Barriers

**Hauptbehauptung:** Übersicht der Ansätze zur RH: Hilbert–Pólya, Berry–Keating xp-Hamilton, GUE-Statistik, Selberg-Spurformel, bekannte Teilresultate (Levinson 1/3, Conrey 2/5), Barrieren.

**Kritische Frage:** RH wird durchgehend als OFFENES PROBLEM (Conjecture) korrekt behandelt. Explizit als Vermutung deklariert, nie als bewiesen. ✅

**Mathematische Korrektheit:** ✅

**Prüfung der Kernaussagen:**

1. **Hilbert–Pólya-Vermutung:** Korrekt als Vermutung formuliert. Beweis-Skizze (selbstadjungiert → reelles Spektrum → Re(ρ)=1/2) logisch korrekt. ✅

2. **Berry–Keating xp-Hamiltonian:**
   - E_n ≈ 2πn/log(n/2πe) stimmt mit γₙ ≈ 2πn/log(n/2π) überein ✅
   - Korrekte Klarstellung: keine rigoros selbstadjungierte Erweiterung gefunden ✅
   - N(E) ≈ E/(2π)·log(E/2πe) + 7/8 mit korrektem Verweis auf Riemann–von Mangoldt ✅

3. **GUE-Definition:** e^(-tr M²)-Dichte ✅. r₂(s) = 1-(sin πs/πs)² ✅

4. **Montgomerys Paarkorrelation:**
   Einschränkung auf Testfunktionen mit f̂ ∈ (-1,1) korrekt vermerkt. ✅
   Bedingte Natur des Satzes (RH-Voraussetzung) explizit. ✅

5. **Selberg-Spurformel:**
   Formel mit λₙ = 1/4 + rₙ², tanh(πr), geodätische Längen korrekt. ✅

6. **Levinson 1/3, Conrey 2/5, Feng 0.41:** Korrekte Attributionen und Jahreszahlen. ✅

7. **Weil-Vermutungen und Deligne:**
   BUG-B5-P24-EN/DE-002 Kommentar im Code weist auf eine frühere Fehlerangabe hin (Weil bewies Vermutungen / Deligne bewies). Die aktuelle Version ist korrekt:
   *"Weil himself proved the special case of curves over finite fields already in 1940; the general case was proved by Deligne (1974)"* ✅
   Kein Fehler mehr vorhanden.

8. **Deligne-Referenz:** \bibitem{Deligne1974} fehlt im Literaturverzeichnis!
   Der Text zitiert `\cite{Deligne1974}` in der Barrier-5-Beschreibung, aber kein entsprechender \bibitem-Eintrag existiert in der bibliography.

**Fehler:**
- [MITTEL] Zeile 356 / Barriere 5: `\cite{Deligne1974}` zitiert, aber kein `\bibitem{Deligne1974}` in der Bibliografie vorhanden. Erzeugt beim LaTeX-Kompilieren eine undefinierte Referenz. → **BUG-B5-P24-EN-001**

**Urteil: ÜBERARBEITUNG NOTWENDIG** (LaTeX-Strukturfehler; Mathematik korrekt)

---

### Paper 24 (DE): Ansätze zur Riemann-Hypothese

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Inhaltlich vollständig übereinstimmend. Alle Formeln und Argumente identisch. ✅

**Fehler:**
- [MITTEL] Zeile 361 / Barriere 5: Analog zu EN fehlt `\bibitem{Deligne1974}`. → **BUG-B5-P24-DE-001**

**Urteil: ÜBERARBEITUNG NOTWENDIG** (LaTeX-Strukturfehler; Mathematik korrekt)

---

## BATCH 6 — ELLIPTISCHE KURVEN / BSD (Papers 25–28)

---

### Paper 25 (EN): Elliptic Curves over Q — Weierstraß Form, Group Law, Torsion, and the Mordell–Weil Theorem

**Hauptbehauptung:** Grundlagen elliptischer Kurven über Q: Weierstraß-Form, j-Invariante, Gruppengesetz (Sekante-Tangente), Mazurs Torsionssatz (15 Typen), Mordell–Weil-Satz, Néron–Tate-Höhe, Reduktion modulo p.

**Kritische Frage:** BSD wird als offenes Problem eingeführt (Remark: *"The Birch–Swinnerton-Dyer conjecture (Paper 27) predicts..."*). Rang-Problem korrekt als ungelöst bezeichnet. ✅

**Mathematische Korrektheit:** ✅

**Prüfung der Kernformeln:**

1. **Weierstraß-Form:** y² = x³ + ax + b, Δ = -16(4a³+27b²) ≠ 0 ✅

2. **j-Invariante:** j(E) = 1728·(4a³)/(4a³+27b²) ✅
   (Anmerkung: Dies ist eine gebräuchliche Form; die vollständige Form wäre
   j = -1728·(4a)³/(Δ/(-16)) = 1728·(4a)³/(4a³+27b²)·(-16/Δ)... Die verwendete Form
   j = 1728·4a³/(4a³+27b²) ist korrekt für das short Weierstraß-Modell.) ✅

3. **Gruppengesetz (Formeln):**
   - Sekantensteigung: λ = (y₂-y₁)/(x₂-x₁) für x₁≠x₂ ✅
   - Tangentensteigung: λ = (3x₁²+a)/(2y₁) für P=Q ✅
   - x₃ = λ²-x₁-x₂ ✅
   - y₃ = λ(x₁-x₃)-y₁ ✅
   Der Beweis via Viète-Formeln ist korrekt: Summe der drei x-Wurzeln = λ². ✅

4. **Inverses:** -(x,y) = (x,-y) ✅ (für short Weierstraß)

5. **Mazurs Torsionssatz:** 15 Gruppen (Z/nZ für n=1,...,10,12 und Z/2Z×Z/2nZ für n=1,2,3,4) ✅

6. **Nagell–Lutz:** x₀,y₀∈Z und y₀²|4a³+27b² ✅

7. **Mordell–Weil:** E(Q) ≅ E(Q)_tors ⊕ Z^r ✅

8. **Néron–Tate-Höhe:** ĥ(P) = lim h(2ⁿP)/4ⁿ ✅
   Bilinearform ⟨P,Q⟩ = ĥ(P+Q)-ĥ(P)-ĥ(Q) ✅

9. **Néron–Ogg–Shafarevich:** E hat gute Reduktion bei p ⟺ Tate-Modul T_ℓ(E) ist bei p unverzweigt. ✅

10. **EN-Beispiel j=1728:** y²=x³-x, Torsion Z/2Z×Z/2Z. Drei 2-Torsionspunkte (0,0),(1,0),(-1,0). ✅
    **HINWEIS:** Der Text sagt "torsion Z/2Z×Z/2Z and an order-4 point over Q(i)" — das ist eine zusätzliche Bemerkung über Q(i), die korrekt ist (4-Torsionspunkt (1,i) über Q(i)), aber nicht mit dem Torsionssatz über Q verwechselt werden sollte. ✅

**Fehler:** Keine mathematischen Fehler.

**Urteil: DRUCKREIF**

---

### Paper 25 (DE): Elliptische Kurven über Q

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Vollständige inhaltliche Übereinstimmung. ✅

**Zusatz in DE:** Beispiel j=0 enthält explizite Nachrechnung j(y²=x³+1) = 1728·0/(0+27) = 0 ✅. Auch extra Angabe der Automorphismen-Ordnung 6 bzw. 4. Dies ist eine sinnvolle Ergänzung.

**Kleiner Unterschied:** DE-Version des Néron–Ogg–Shafarevich-Satzes enthält die Galois-Darstellung explizit (ρ_ℓ: Gal(Q̄/Q) → GL₂(Z_ℓ)), EN-Version etwas kompakter. Beide korrekt. ✅

**Fehler:** Keine.

**Urteil: DRUCKREIF**

---

### Paper 26 (EN): The L-Function of an Elliptic Curve — Euler Product, Frobenius Trace, Analytic Continuation, and the Functional Equation

**Hauptbehauptung:** L-Funktion L(E,s): Frobenius-Spur aₚ, Hasse-Schranke, Euler-Produkt, Dirichlet-Reihe, vervollständigte Funktion, Funktionalgleichung via Modularitätssatz, Wiles/BCDT.

**Kritische Frage:** BSD wird nur als Vermutung zitiert: *"The BSD conjecture predicts r_an = r"* — BSD korrekt als offen behandelt. Modularitätssatz korrekt als bewiesen dargestellt (Wiles 1995 / BCDT 2001). ✅

**Mathematische Korrektheit:** ✅

**Prüfung der Kernformeln:**

1. **Frobenius-Spur:** aₚ = p+1-#E(𝔽ₚ) ✅

2. **Hasse-Schranke:** |aₚ| ≤ 2√p ✅
   Beweis via α·ᾱ = p korrekt. ✅

3. **Lokale L-Faktoren:**
   - Gute Reduktion: L_p(E,s)⁻¹ = 1-aₚp⁻ˢ+p^(1-2s) ✅
   - Multiplikative Reduktion: L_p(E,s)⁻¹ = 1-ε_p p⁻ˢ, ε_p∈{±1} ✅
   - Additive Reduktion: L_p(E,s)⁻¹ = 1 ✅

4. **Euler-Produkt:** Konvergenz für Re(s)>3/2 durch |aₙ|≤(n+1)p^(n/2) ✅

5. **Dirichlet-Reihen-Rekursion:** a_{p^k} = aₚ·a_{p^(k-1)} - p·a_{p^(k-2)} ✅

6. **Vervollständigte L-Funktion:**
   Λ(E,s) = (√N_E/(2π))^s · Γ(s) · L(E,s) ✅

7. **Funktionalgleichung:** Λ(E,s) = ε(E)·Λ(E,2-s), ε(E)∈{±1} ✅

8. **Wurzelzahl und Parität:**
   ε(E)=-1 → Λ(E,1) = -Λ(E,1) → L(E,1)=0 ✅
   Paritätsvermutung (-1)^r = ε(E) korrekt als Vermutung deklariert. ✅

9. **Beispiel y²=x³-x, p=5:**
   Punkte: (0,0),(1,0),(2,1),(2,4),(3,2),(3,3),(4,0) → 7 affin + O = 8
   a₅ = 5+1-8 = -2, |-2| ≤ 2√5 ≈ 4.47 ✅

10. **Beispiel y²=x³-25x, Rang 1:**
    Punkt (-4,6): 6²=36, (-4)³-25(-4)=-64+100=36 ✅

**Fehler:** Keine mathematischen Fehler.

**Urteil: DRUCKREIF**

---

### Paper 26 (DE): Die L-Funktion einer elliptischen Kurve

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Vollständige Übereinstimmung aller Formeln. ✅

**Zusatz in DE:** Numerischer Abschnitt (Sec. 7) ausführlicher als EN:
- Rang 0: E₀: y²=x³-x, N=32, Frobenius-Spuren a₃=0,a₅=-2,a₇=0,... explizit
- L(E₀,1) ≈ 0.6551, Kolyvagin-Beweis erwähnt ✅
- Rang 1: E₁: y²+y=x³-x (N=37), a₂=-2,a₃=-3,..., L'(E₁,1)≈0.3059, ε(E₁)=-1 ✅
  **WICHTIG:** Diese Kurve y²+y=x³-x ist nicht in short Weierstraß-Form; der DE-Text behandelt sie korrekt ohne die short-form-Beschränkung. ✅

**Fehler:** Keine.

**Urteil: DRUCKREIF**

---

### Paper 27 (EN): The Birch–Swinnerton-Dyer Conjecture

**Hauptbehauptung:** Formulierung der schwachen und starken BSD-Vermutung, Tate–Shafarevich-Gruppe, bekannte Fälle (Coates–Wiles, Kolyvagin), numerische Evidenz, Abstieg/Selmer-Gruppe.

**Kritische Frage:** BSD wird durchgehend als OFFENES PROBLEM (Conjecture) korrekt behandelt. Explizit als Millennium-Preisaufgabe bezeichnet. Schwache BSD nur für analytischen Rang ≤ 1 bewiesen. ✅

**Mathematische Korrektheit:** ✅

**Prüfung der Kernformeln:**

1. **Schwache BSD:** ord_{s=1} L(E,s) = r ✅

2. **Starke BSD-Formel:**
   `lim_{s→1} L(E,s)/(s-1)^r = (Ω_E · |Sha(E/Q)| · R_E · ∏_p cₚ) / |E(Q)_tors|²` ✅
   Alle Zutaten korrekt definiert (Période Ω_E, Néron–Tate-Regulator R_E, Tamagawa-Zahlen cₚ). ✅

3. **BSD-Heuristikprodukt:** B_E(X) = ∏_{p≤X, good} #E(𝔽ₚ)/p ~ C_E(log X)^r ✅

4. **Tate–Shafarevich-Gruppe:**
   Sha(E/Q) = ker(H¹(Q,E) → ∏_v H¹(Q_v,E)) ✅
   Cassels–Tate-Paarung alternierend → |Sha| ist Quadrat ✅

5. **Coates–Wiles (1977):** Korrekt zitiert für CM-Kurven. ✅

6. **Kolyvagin (1990):**
   (i) L(E,1)≠0 → r=0 und Sha endlich ✅
   (ii) L'(E,1)≠0 → r=1 und Sha endlich ✅
   Beweis-Skizze via Euler-Systeme von Heegner-Punkten und Gross–Zagier korrekt. ✅

7. **Exakte Folge für Selmer-Gruppe:**
   0 → E(Q)/2E(Q) → Sel₂(E/Q) → Sha(E/Q)[2] → 0 ✅

8. **Dimensionsformel:**
   `r + dim_{F₂} Sha[2] = dim_{F₂} Sel₂(E/Q) - dim_{F₂} E(Q)_tors[2]` ✅

**Fehler:** Keine mathematischen Fehler.

**Urteil: DRUCKREIF**

---

### Paper 27 (DE): Die Birch–Swinnerton-Dyer-Vermutung

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Vollständige Übereinstimmung. ✅

**Zusatz in DE:** Abschnitt 6 (Abstieg und Selmer-Gruppe) in DE deutlich ausführlicher: Kummer-Einbettung erklärt, strukturelle Herleitung der Dimensionsformel aus der exakten Folge und dem Struktursatz E(Q) ≅ Z^r × E(Q)_tors explizit. Diese Ergänzungen sind mathematisch korrekt und didaktisch wertvoll. ✅

**Fehler:** Keine.

**Urteil: DRUCKREIF**

---

### Paper 28 (EN): Congruent Numbers, Elliptic Curves, and BSD

**Hauptbehauptung:** n kongruent ⟺ Rang(E_n) ≥ 1, mit E_n: y²=x³-n²x. Torsion = Z/2Z×Z/2Z. Tunnells Satz (bedingt unter BSD). Verifikation für n=5,6,7. Status n=1,2,3.

**Kritische Frage:** BSD wird als Annahme für Tunnells Umkehrung klar deklariert: *"Conversely (assuming BSD for E_n)"*. ✅ Das kongruente Zahlenproblem korrekt als "conditionally resolved" beschrieben. ✅

**Mathematische Korrektheit:** ✅

**Prüfung der Kernformeln und -beweise:**

1. **Äquivalenz-Beweis (Theorem 2.1):**
   Parametrisierung x₀ = c²/4, y₀ = c(a²-b²)/8.
   Verifikation: c⁴-16n² = (a²+b²)²-(ab)² = (a²-b²)² ✅
   Daher: x₀³-n²x₀ = c²(a²-b²)²/64 = y₀² ✅

   Rückrichtung: a = |（x₀²-n²)/y₀|, b = |2nx₀/y₀|, c = |(x₀²+n²)/y₀|
   ab/2 = nx₀(x₀²-n²)/y₀² = n (da y₀²=x₀(x₀²-n²)) ✅

   **HINWEIS:** Der Schritt x₀(x₀²-n²) = x₀³-n²x₀ = y₀² setzt voraus, dass x₀≠0.
   Für x₀=0 wäre y₀=0 (Torsionspunkt). Der Beweis nimmt implizit y₀≠0 an, was durch
   "point of infinite order" gesichert ist. ✅

2. **Torsion von E_n:**
   E_n: y²=x(x-n)(x+n) hat 2-Torsionspunkte (0,0),(n,0),(-n,0) ✅
   Keine weiteren Torsionspunkte via Nagell–Lutz: y₀²|4a³+27b² = 4(-n²)³+0 = -4n⁶, also y₀²|4n⁶ ✅

3. **Tunnells Satz:**
   Formeln für A(n), B(n), C(n), D(n) korrekt aufgeführt. ✅
   Beweisrahmen (Modulformen Gewicht 3/2, Shimura-Lifting, Theta-Funktionen) korrekt beschrieben. ✅
   Die konditionelle Natur der Umkehrung klar. ✅

4. **Verifikation n=6:** Punkt (25/4, 35/8) auf E₆: y²=x³-36x
   (35/8)²=1225/64; (25/4)³-36(25/4)=15625/64-14400/64=1225/64 ✅

5. **Verifikation n=5:** Dreieck (3/2, 20/3, 41/6), Fläche=5.
   (3/2)²+(20/3)²=9/4+400/9=81/36+1600/36=1681/36=(41/6)² ✅
   Fläche=(1/2)(3/2)(20/3)=5 ✅
   Punkt (25/4, 75/8) auf E₅: (75/8)²=5625/64; (25/4)³-25(25/4)=15625/64-10000/64=5625/64 ✅

6. **Verifikation n=7:** (35/12)²+(24/5)²=1225/144+576/25=(30625+82944)/3600=113569/3600=(337/60)² ✅
   Fläche=(1/2)(35/12)(24/5)=(1/2)(840/60)=(1/2)(14)=7 ✅

7. **Status n=1,2,3:**
   - n=1: A(1)=2, B(1)=2 → A(1)=2 ≠ 4=2B(1) → nicht kongruent unter BSD.
     **ANMERKUNG:** Der Text schreibt "A(1) = 2, B(1) = 2, so A(1)=2 ≠ 4=2B(1)". Das sollte
     gelesen werden als: A(1)=2 und 2B(1)=4, also A(1)≠2B(1). Korrekt. ✅
     Fermat-Beweis (unendlicher Abstieg) wird als unabhängiger unbedingter Beweis erwähnt. ✅
   - n=2,3: Korrekt als "not yet unconditionally proven" deklariert. ✅

**Fehler:** Keine mathematischen Fehler.

**Urteil: DRUCKREIF**

---

### Paper 28 (DE): Kongruente Zahlen, elliptische Kurven und BSD

**Mathematische Korrektheit:** ✅

**EN/DE-Parität:** Vollständige inhaltliche Übereinstimmung. Alle Rechnungen identisch. ✅

**Kleine strukturelle Abweichung:** DE-Version referenziert Kobyvagin im Status n=1 nicht explizit. Die EN-Version nennt Fermat auch. Kein inhaltlicher Fehler, nur Stilunterschied.

**Fehler:** Keine.

**Urteil: DRUCKREIF**

---

## ZUSAMMENFASSUNG ALLER NEUEN BUGS

| Bug-ID | Schwere | Paper | Zeile | Beschreibung |
|---|---|---|---|---|
| BUG-B5-P24-EN-001 | MITTEL | Paper 24 EN | Z. 356 | \cite{Deligne1974} zitiert, aber kein \bibitem vorhanden |
| BUG-B5-P24-DE-001 | MITTEL | Paper 24 DE | Z. 361 | Gleicher fehlender \bibitem{Deligne1974} |

---

## GESAMTÜBERSICHT BATCH 5 & 6

| Paper | Titel | EN | DE | Urteil |
|---|---|---|---|---|
| 21 | Riemann Zeta-Funktion | ✅ DRUCKREIF | ✅ DRUCKREIF | Mathematisch makellos |
| 22 | Nichttriviale Nullstellen & GUE | ✅ DRUCKREIF | ✅ DRUCKREIF | Mathematisch makellos |
| 23 | Explizite Formel & PZS | ✅ DRUCKREIF | ✅ DRUCKREIF | s=0-Beitrag korrekt erklärt |
| 24 | Ansätze zur RH | ⚠️ ÜBERARBEITUNG | ⚠️ ÜBERARBEITUNG | Fehlende Deligne-Referenz |
| 25 | Elliptische Kurven über Q | ✅ DRUCKREIF | ✅ DRUCKREIF | Mathematisch makellos |
| 26 | L-Funktion ell. Kurven | ✅ DRUCKREIF | ✅ DRUCKREIF | Modularitätssatz korrekt |
| 27 | BSD-Vermutung | ✅ DRUCKREIF | ✅ DRUCKREIF | Vollständige Darstellung |
| 28 | Kongruente Zahlen & BSD | ✅ DRUCKREIF | ✅ DRUCKREIF | Alle Rechnungen verifiziert |

### Kritische Sicherheitsfragen (RH/BSD als bewiesen?)

- **RH** wird in ALLEN 8 Riemann-Papers korrekt als offenes Problem/Vermutung behandelt. Kein Paper behauptet einen Beweis. ✅ KEINE FATALEN FEHLER.
- **BSD** wird in ALLEN 8 BSD/Elliptische-Kurven-Papers korrekt als Vermutung behandelt. Weder schwache noch starke BSD wird als bewiesen deklariert (nur die partiellen Resultate von Coates–Wiles und Kolyvagin für Rang ≤ 1 werden korrekt als bewiesen dargestellt). ✅ KEINE FATALEN FEHLER.

### Mathematische Qualität

Die Papers zeigen insgesamt ein sehr hohes mathematisches Niveau:
- Alle Funktionalgleichungen korrekt
- Alle Beweis-Skizzen mathematisch valide
- Alle Formeln verifiziert (ψ-Formel, BSD-Formel, Gruppengesetz, Tunnell-Formeln)
- Zahlenwerte numerisch korrekt (erste Nullstellen, a₅-Berechnung, Verifikationen n=5,6,7)
- EN/DE-Parität in allen Fällen gewahrt

### Einzige Mängel

Paper 24 (EN+DE): fehlende Deligne-Referenz — leicht behebbar, mathematisch unbedenklich.
