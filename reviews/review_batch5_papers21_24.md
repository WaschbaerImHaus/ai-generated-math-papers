# Mathematisches Gutachten — Batch 5: Papers 21–24 (Riemann-Hypothese)

**Datum:** 2026-03-11
**Gutachter:** Claude Sonnet 4.6
**Geprüfte Dateien:**
- `papers/batch5/paper21_riemann_zeta_function_en.tex`
- `papers/batch5/paper21_riemann_zeta_function_de.tex`
- `papers/batch5/paper22_nontrivial_zeros_en.tex`
- `papers/batch5/paper22_nontrivial_zeros_de.tex`
- `papers/batch5/paper23_explicit_formula_pnt_en.tex`
- `papers/batch5/paper23_explicit_formula_pnt_de.tex`
- `papers/batch5/paper24_rh_approaches_en.tex`
- `papers/batch5/paper24_rh_approaches_de.tex`

---

## Gesamturteil

Alle vier Paper-Paare (EN+DE) sind inhaltlich sehr solide und entsprechen dem Stand
der Fachliteratur. Die mathematischen Kernaussagen (Dirichlet-Reihe, Euler-Produkt,
Funktionalgleichung, Riemann-von-Mangoldt-Formel, explizite Formel, RH-Ansätze)
sind korrekt formuliert. Es gibt mehrere Bugs unterschiedlicher Schwere.

---

## Paper 21: Die Riemann'sche Zeta-Funktion (EN + DE)

### Mathematische Korrektheit

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| Absolute Konvergenz für Re(s)>1 | ✓ Korrekt | Beweis via Weierstraß M-Test vollständig |
| Euler-Produkt | ✓ Korrekt | Beweis via FTA klar |
| Gamma-Funktion Eigenschaften | ✓ Korrekt | Alle 4 Eigenschaften korrekt |
| Γ(s)ζ(s) = ∫t^{s-1}/(e^t-1)dt | ✓ Korrekt | Beweis via dominierte Konvergenz korrekt |
| Pol bei s=1, Residuum 1 | ✓ Korrekt | Laurent-Entwicklung korrekt |
| Funktionalgleichung (eq:functional) | ✓ Korrekt | Via Jacobi-Theta, Beweisskizze ausreichend |
| ξ(s) = ξ(1-s) | ✓ Korrekt | |
| Triviale Nullstellen bei s=-2,-4,... | ✓ Korrekt | Beweis via sin(-mπ)=0 vollständig |
| Riemann-Siegel-Formel Z(t) | ✓ Korrekt | |
| Euler-Formel ζ(2k) | ✓ Korrekt | Formel mit (-1)^{k+1} korrekt für k≥1 |

### LaTeX-Korrektheit

**EN-Version:**
- Alle `\newtheorem`-Umgebungen korrekt definiert
- `\Q` fehlt im Preamble (aber `\Q` wird im Text nicht verwendet — kein Fehler)
- Alle `\ref` verweisen auf tatsächlich definierte Labels
- `\bibitem` für alle `\cite` vorhanden: Riemann1859, Davenport2000, Titchmarsh1986,
  Ahlfors1979, Apostol1976, Platt2021 ✓

**DE-Version:**
- `\usepackage[ngerman]{babel}` vorhanden ✓
- `\Q` im Preamble definiert ✓ (wird in Abschnitt 4 nicht direkt gebraucht)
- Alle `\bibitem` vorhanden ✓

### EN↔DE Konsistenz

- Abschnittsstruktur identisch ✓
- Formelbezeichnungen identisch (`eq:functional`) ✓
- **BUG-B5-01:** DE-Version enthält in der θ(t)-Asymptotik (Def. der Riemann-Siegel-Theta)
  nur `+ O(t^{-3})`, während die EN-Version das präzisere `- 7/(5760 t^3) + O(t^{-5})` angibt.
  Das DE-Paper unterschlägt einen Term und ist weniger präzise.

### Gefundene Bugs

**BUG-B5-01** (MINOR)
**Datei:** `paper21_riemann_zeta_function_de.tex`, Zeile 362
**Problem:** Asymptotik der Riemann-Siegel-Theta-Funktion im DE-Paper unvollständig.
EN-Version: `- \frac{7}{5760\,t^3} + O(t^{-5})`; DE-Version: `+ O(t^{-3})`.
Dies ist inhaltlich nicht falsch (korrekte O-Ordnung), aber inkonsistent mit der EN-Version.

---

## Paper 22: Die nicht-trivialen Nullstellen (EN + DE)

### Mathematische Korrektheit

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| Riemann-von-Mangoldt-Formel N(T) | ✓ Korrekt | Konstante 7/8 korrekt |
| Erste Nullstelle γ₁ = 14.134725... | ✓ Korrekt | Wert korrekt |
| Die ersten 10 Nullstellen (γ₁–γ₁₀) | ✓ Korrekt | Alle Werte auf 6 Stellen korrekt |
| Gram-Punkte: Erstversagen bei n=126 | ✓ Korrekt | Historisch korrekt (Hutchinson 1925, Lehmer 1956) |
| Normalisierter Abstand δ_n | ✓ Korrekt | Formel korrekt |
| Montgomery-Paarkorrelation | ✓ Korrekt | GUE-Funktion 1-(sin πx/πx)² korrekt |
| 10^{20}-te Nullstelle ≈ 1.52×10^{19} | ✓ Korrekt | Odlyzko 1987 korrekt |
| Platt-Trudgian 2021: T=3×10^{12} | ✓ Korrekt | |

### LaTeX-Korrektheit

**EN-Version:**
- `\Q` fehlt in Custom commands (wird im Text nicht benötigt — kein Fehler)
- Alle `\bibitem` vorhanden ✓
- `\url` für LMFDB: `hyperref` ist geladen ✓

**BUG-B5-02** (MITTEL)
**Datei:** `paper22_nontrivial_zeros_en.tex`, Zeile 69
**Problem:** Im Einleitungsabsatz steht: `the functional equation ζ(s) = ζ(1-s) (up to explicit factors)`.
Das ist ungenau bis missverständlich. Die korrekte Gleichung ist
`ζ(s) = 2^s π^{s-1} sin(πs/2) Γ(1-s) ζ(1-s)`.
Die Formulierung „up to explicit factors" verschleiert, dass die Faktoren wesentlich sind.
In Paper 21 ist die Gleichung korrekt formuliert. Die DE-Version formuliert dies besser
(„Aus der Funktionalgleichung folgt: Ist ρ eine Nullstelle, so auch 1-ρ") ohne die
fehlerhafte Vereinfachung.

**DE-Version:**
- Diese missverständliche Vereinfachung taucht in der DE-Version nicht auf ✓
- `\bibitem{LMFDB}` vorhanden ✓

### EN↔DE Konsistenz

- Alle Numerics (γ₁–γ₁₀) konsistent zwischen EN und DE ✓
- DE verwendet europäisches Dezimalkomma (14,134725...) — korrekt für DE-Dokument ✓
- Tabellenformat vs. align*-Umgebung für Nullstellen: beide Versionen verwenden align* ✓

---

## Paper 23: Explizite Formel und Primzahlsatz (EN + DE)

### Mathematische Korrektheit

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| von-Mangoldt-Funktion Λ(n) | ✓ Korrekt | |
| Tschebyschew ψ(x) = θ(x) + θ(x^{1/2}) + ... | ✓ Korrekt | |
| π(x) = θ(x)/log x + ∫θ(t)/(t(log t)²)dt | ✓ Korrekt | Abel-Summation korrekt |
| -ζ'(s)/ζ(s) = Σ Λ(n)/n^s | ✓ Korrekt | Beweis via Euler-Produkt vollständig |
| Perrons Formel | ✓ Korrekt | |
| Explizite Formel für ψ(x) | ✓ Korrekt | ζ'(0)/ζ(0) = log(2π) korrekt |
| Primzahlsatz via 3+4cosθ+cos2θ≥0 | ✓ Korrekt | Mertens-Argument korrekt |
| Fehlerterm O(x^{β₀}(log x)²) | ✓ Korrekt | |
| Schoenfeld-Schranke unter RH | ✓ Korrekt | Konstante 1/(8π) korrekt |
| Riemanns explizite Formel für π(x) | ✓ Korrekt | J(x)-Ansatz korrekt |

**BUG-B5-03** (MINOR)
**Datei:** `paper23_explicit_formula_pnt_en.tex`, Zeile 318–319
**Problem:** In der Bemerkung nach Theorem 7.1 steht:
`The integral in (eq:riemann_explicit) equals O(x^{-1} log x)^{-1} as x → ∞`.
Das ist ein Tippfehler/Klammerungsfehler: `O(x^{-1}\log x)^{-1}` sollte
`O((x \log x)^{-1})` oder schlicht `O(1/\log x)` sein (was auch im Abstract steht).
Die -1-Potenz steht außerhalb der O-Notation, das ist syntaktisch falsch.

### LaTeX-Korrektheit

**EN + DE:** Alle Umgebungen korrekt, alle `\bibitem` vorhanden ✓
DE enthält `\DeclareMathOperator{\Li}{Li}` (wird in Formeln verwendet) ✓

### EN↔DE Konsistenz

- Beide Versionen konsistent ✓
- Schoenfeld-Wert 2657 in beiden Versionen korrekt ✓
- Beweisskizze PZS: EN-Version erwähnt zusätzlich `Dirichlet's theorem with error term under GRH`
  als eigenständigen Abschnitt — die DE-Version ist kürzer, aber nicht inkonsistent.
- **BUG-B5-04** (MINOR): Der Remark im DE-Paper nach Thm. für ψ(x)-Fehlerterm
  enthält einen Schreibfehler: `die gesamt Summe` statt `die gesamte Summe` (Zeile ~263 DE).

---

## Paper 24: Ansätze zur Riemann-Hypothese (EN + DE)

### Mathematische Korrektheit

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| Hilbert-Pólya: Eigenwerte ⟹ RH | ✓ Korrekt | Beweis trivial aber korrekt |
| Berry-Keating H=xp | ✓ Korrekt | |
| E_n ≈ 2πn/log(n/2πe) stimmt mit γ_n | ✓ Korrekt | Asymptotik konsistent mit Riemann-von-Mangoldt |
| N(E) ≈ E/(2π)·log(E/2πe)+7/8 | ✓ Korrekt | Konsistent mit Paper 22 |
| Montgomery-Paarkorrelation (GUE) | ✓ Korrekt | Vorbehalt (Fourier-Träger) korrekt erwähnt |
| Selberg-Spurformel | ✓ Korrekt | |
| Levinson 1/3, Conrey 2/5, Feng 0.41 | ✓ Korrekt | |
| De la Vallée-Poussin: nullstellenfreies Gebiet σ > 1 - c/log(|t|+2) | ✓ Korrekt | |
| Weil (1940) + Deligne (1974): Riemann-Hyp. für endl. Körper | ✓ Korrekt | |

**BUG-B5-05** (MINOR)
**Datei:** `paper24_rh_approaches_de.tex`, Zeile 184
**Problem:** Im Diskussionsteil des Berry-Keating-Theorems steht:
`keiner liefert bisher einen vollst\"andigen rigoros Beweis`.
Grammatikfehler: `rigoros` muss `rigorosen` heißen (Adjektiv-Deklination).

**BUG-B5-06** (MITTEL)
**Datei:** `paper24_rh_approaches_de.tex`, Zeile 204 (GUE-Definition)
**Problem:** In der GUE-Definition wird `\tr` als `\DeclareMathOperator{\tr}{Spur}` definiert,
aber im Text ist die Dichte proportional zu `e^{-\tr M^2}`. Die Spur im mathematischen Kontext
heißt standardmäßig `tr` (nicht `Spur`) — auch für ein deutsches Dokument, da es sich um
mathematischen Fachbegriff handelt. Dies führt beim Kompilieren zu `e^{-\text{Spur}(M^2)}`
statt `e^{-\mathrm{tr}(M^2)}`. In der EN-Version ist `\tr` korrekt als `tr` definiert.

### LaTeX-Korrektheit

**EN-Version:**
- `\conjecture` Umgebung korrekt definiert ✓
- `\spec` als Operator korrekt deklariert ✓
- Alle `\bibitem` vorhanden: Riemann1859, Montgomery1973, Odlyzko1987, BerryKeating1999,
  Connes1999, Conrey1989, Davenport2000, Titchmarsh1986 ✓

**DE-Version:**
- `\conjecture` als `Vermutung` korrekt ✓
- `\spec` als `Spec` — abweichend von EN (`spec`); funktioniert aber ✓
- `\tr` als `Spur` (BUG-B5-06)

### EN↔DE Konsistenz

- Inhalt vollständig äquivalent ✓
- Barrieren-Abschnitt (5 Barrieren) konsistent ✓
- 3-Wege-Analogie-Tabelle in beiden Versionen vorhanden ✓

---

## Zusammenfassung Batch 5

| Bug-ID | Schwere | Datei | Beschreibung |
|--------|---------|-------|--------------|
| BUG-B5-01 | MINOR | paper21_de.tex | θ(t)-Asymptotik im DE-Paper unvollständig (Terme nach O(t^{-3}) fehlen) |
| BUG-B5-02 | MITTEL | paper22_en.tex | Funktionalgleichung als „ζ(s)=ζ(1-s) up to factors" vereinfacht |
| BUG-B5-03 | MINOR | paper23_en.tex | Klammerfehler in O-Notation: `O(x^{-1}\log x)^{-1}` statt `O((x\log x)^{-1})` |
| BUG-B5-04 | MINOR | paper23_de.tex | Grammatik: „die gesamt Summe" statt „die gesamte Summe" |
| BUG-B5-05 | MINOR | paper24_de.tex | Grammatik: „rigoros Beweis" statt „rigorosen Beweis" |
| BUG-B5-06 | MITTEL | paper24_de.tex | `\DeclareMathOperator{\tr}{Spur}` produziert Spur statt tr in Formeln |

**Keine fundamentalen mathematischen Fehler gefunden.** Alle Hauptsätze sind korrekt,
vollständig und gut referenziert. Die LaTeX-Struktur ist durchgehend kompilierbar.
