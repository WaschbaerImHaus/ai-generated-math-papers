# Mathematisches Gutachten — Batch 6: Papers 25–28 (Elliptische Kurven + BSD)

**Datum:** 2026-03-11
**Gutachter:** Claude Sonnet 4.6
**Geprüfte Dateien:**
- `papers/batch6/paper25_elliptic_curves_Q_en.tex`
- `papers/batch6/paper25_elliptic_curves_Q_de.tex`
- `papers/batch6/paper26_l_function_elliptic_en.tex`
- `papers/batch6/paper26_l_function_elliptic_de.tex`
- `papers/batch6/paper27_bsd_conjecture_en.tex`
- `papers/batch6/paper27_bsd_conjecture_de.tex`
- `papers/batch6/paper28_congruent_numbers_bsd_en.tex`
- `papers/batch6/paper28_congruent_numbers_bsd_de.tex`

---

## Gesamturteil

Papers 25–27 sind mathematisch korrekt und vollständig. Paper 28 (kongruente Zahlen)
enthält einen schwerwiegenden inhaltlichen Fehler im Beweis (unvollständige/abgebrochene
Beweisführung mit sichtbaren Entwurfsresten) sowie einen numerischen Fehler in der
Diskussion der Tunnell-Kriterien für kleine n.

---

## Paper 25: Elliptische Kurven über Q (EN + DE)

### Mathematische Korrektheit

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| Weierstraß-Form, Diskriminante Δ=-16(4a³+27b²) | ✓ Korrekt | |
| j-Invariante = 1728·4a³/(4a³+27b²) | ✓ Korrekt | |
| Isomorphieklassifikation durch j | ✓ Korrekt | „bis auf Twist" über Q korrekt erwähnt |
| Gruppengesetz (Sekanten-Tangenten) | ✓ Korrekt | Formeln λ, x₃, y₃ korrekt |
| Viète-Argument für x₃=λ²-x₁-x₂ | ✓ Korrekt | |
| Mazur-Torsionssatz: 15 Typen | ✓ Korrekt | |
| Nagell-Lutz-Satz | ✓ Korrekt | |
| Mordell-Weil: E(Q) ≅ E(Q)_tors ⊕ Z^r | ✓ Korrekt | |
| Schwacher Mordell-Weil | ✓ Korrekt | 2-Abstieg korrekt skizziert |
| Kanonische Höhe N̂: ĥ(nP) = n²ĥ(P) | ✓ Korrekt | Grenzwert-Def. korrekt |
| Gute/schlechte Reduktion | ✓ Korrekt | |
| Néron-Ogg-Shafarevich | ✓ Korrekt | Statement korrekt |

**BUG-B6-01** (MINOR)
**Datei:** `paper25_elliptic_curves_Q_en.tex`, Zeilen 100–103
**Problem:** Die j-Invariante wird doppelt, inkonsistent angegeben:
```
j(E) = -1728 · (4a)^3/Δ  (Zeile 100-101, erste Formel)
j(E) = 1728 · 4a³/(4a³+27b²)  (Zeile 101, nach `=`)
```
Die erste Fassung `j(E) = -1728 · (4a)^3/Δ` mit Δ = -16(4a³+27b²) ist korrekt
(da Δ negativ-vorzeichenbehaftet), aber die Zwischendarstellung `= -1728 · 4a³/(Δ/(-16))·(-16)`
ist verwirrend und enthält redundante Rechenschritte. In der DE-Version steht nur
`j(E) = 1728 · 4a³/(4a³+27b²)` — das ist die sauberste Form.
Die EN-Version sollte auf die einfache Form vereinheitlicht werden.

**BUG-B6-02** (MINOR)
**Datei:** `paper25_elliptic_curves_Q_en.tex`, Zeile 233–234
**Problem:** Im Example-Block zum singulären Fall steht:
`E: y^2 = x^3 - x^2 (with b=0, so Δ=0) is singular`.
Das ist korrekt (`b=0` und `a=-1, b=0` gibt Δ=-16(4(-1)³+27·0²)=-16·(-4)=64 ≠ 0),
aber hier ist `a=0, b=0` gemeint (für y²=x³-x² müsste man schreiben
y²=x³+0·x+0 nach Umformung, die aber nicht direkt in Weierstraß-Kurzform ist).
Tatsächlich ist y²=x³-x² eine singuläre Kurve (Knoten bei x=0) — aber sie hat
`a=0, b=0` in Weierstraß-Kurzform (nach x→x+1/3-Verschiebung ist sie y²=x³-1/3·x-2/27).
Die Angabe `b=0` ist falsch: die Diskriminante der Form y²=x³-x² ist Δ=0 wegen
der Singularität, aber nicht weil b=0. Dies ist ein Darstellungsfehler.

### LaTeX-Korrektheit

**EN-Version:**
- `\gcd` als Operator deklariert — aber `\gcd` ist schon in LaTeX als `\gcd` vordefiniert;
  `\DeclareMathOperator{\gcd}{gcd}` ist redundant, aber nicht fehlerhaft ✓
- `\Fp` korrekt definiert als `\mathbb{F}_p` ✓
- Alle `\bibitem` vorhanden: Silverman1986, Cassels1991, Mazur1977, Cornell1997, Washington2008 ✓
- `\example` Umgebung ist NICHT in der EN-Version deklariert, wird aber als `\begin{example}` verwendet!

**BUG-B6-03** (HOCH)
**Datei:** `paper25_elliptic_curves_Q_en.tex`
**Problem:** Die `example`-Umgebung wird zweimal verwendet (`\begin{example}` in Zeile 123
und Zeile 229), ist aber im Preamble NICHT deklariert. Die DE-Version hat:
`\newtheorem{example}[theorem]{Beispiel}`, die EN-Version fehlt:
`\newtheorem{example}[theorem]{Example}`.
Dies führt zu einem LaTeX-Kompilierfehler (Undefined environment `example`).

**DE-Version:**
- `\newtheorem{example}` korrekt definiert ✓
- `\bibitem{Cornell1997}` fehlt in der DE-Version!

**BUG-B6-04** (MITTEL)
**Datei:** `paper25_elliptic_curves_Q_de.tex`
**Problem:** `\bibitem{Cornell1997}` fehlt in der Bibliographie der DE-Version,
obwohl es in der EN-Version vorhanden ist. Da keine Referenz `\cite{Cornell1997}`
im DE-Text vorkommt (der DE-Text referenziert Cornell1997 nicht), handelt es sich
technisch um kein LaTeX-Fehler — aber die Bibliographien sind inkonsistent.

### EN↔DE Konsistenz

- Mordell-Weil-Beweis: EN hat präzisere Formulierungen (Parallelogramm-Gesetz mit eq-labels) ✓
- DE-Version fehlt das Néron-Ogg-Shafarevich-Theorem als eigener Satz — die DE-Version
  endet nach der Reduktionsdefinition ohne den abschließenden Theorem-Block.

**BUG-B6-05** (MITTEL)
**Datei:** `paper25_elliptic_curves_Q_de.tex`
**Problem:** Das Theorem `Néron-Ogg-Shafarevich` (Theorem 6.2 im EN-Paper) fehlt
vollständig in der DE-Version. Der Reductions-Abschnitt endet nach der Definition,
ohne den Satz zu formulieren.

---

## Paper 26: L-Funktion elliptischer Kurven (EN + DE)

### Mathematische Korrektheit

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| Frobenius-Spur a_p = p+1-#E(F_p) | ✓ Korrekt | |
| Hasse: |a_p| ≤ 2√p | ✓ Korrekt | Beweis via Eigenwerte korrekt |
| Lokale L-Faktoren (gute/mult./add. Red.) | ✓ Korrekt | |
| L(E,s) konvergiert für Re(s)>3/2 | ✓ Korrekt | |
| Rekursion a_{p^k} = a_p a_{p^{k-1}} - p a_{p^{k-2}} | ✓ Korrekt | |
| Vervollständigte L-Funktion Λ(E,s) | ✓ Korrekt | |
| Funktionalgleichung Λ(E,s)=ε(E)Λ(E,2-s) | ✓ Korrekt | Via Modularitätssatz |
| Modulärer Isomorphismus X₀(N)→E | ✓ Korrekt | |
| Fermats Letzter Satz als Korollar | ✓ Korrekt | Via Ribet-Satz korrekt erwähnt |
| Wurzelzahl ε(E)=-1 ⟹ L(E,1)=0 | ✓ Korrekt | |

**Numerische Prüfung Beispiel a₅ für y²=x³-x:**
- x=0: x³-x=0 ✓, y=0 ✓ (1 Punkt)
- x=1: 1-1=0 ✓, y=0 ✓ (1 Punkt)
- x=2: 8-2=6≡1 (mod 5) ✓, y²=1 → y=1,4 ✓ (2 Punkte)
- x=3: 27-3=24≡4 (mod 5) ✓, y²=4 → y=2,3 ✓ (2 Punkte)
- x=4: 64-4=60≡0 (mod 5) ✓, y=0 ✓ (1 Punkt)
Gesamt: 7 affine + O = 8, a₅=5+1-8=-2 ✓

**BUG-B6-06** (MITTEL)
**Datei:** `paper26_l_function_elliptic_en.tex`, Zeilen 308–312 (Beispiel 6.2)
**Problem:** Im Beispiel für y²=x³-25x, rank-1-Kurve, wird zunächst behauptet,
der Generator unendlicher Ordnung sei „(5,0)", dann korrigiert: „(5,0) ist ein
2-Torsionspunkt". Die eigentliche Korrektur gibt den Generator als `(-4, 6)` an und
verifiziert: 6²=36, (-4)³-25(-4)=-64+100=36 ✓. Das ist korrekt.
JEDOCH: Der Generator der freien Teils von E(Q) für y²=x³-25x ist tatsächlich
(-4, 6)? Prüfung: y²=x³-25x bei x=-4: (-4)³-25·(-4)=-64+100=36=6². ✓
Der Text enthält den Entwurfsüberrest „(5, 0)... wait," — das sind sichtbare
Bearbeitungsspuren, die aus dem Endprodukt entfernt werden müssen.

### LaTeX-Korrektheit

**EN + DE:**
- Alle Umgebungen korrekt definiert ✓
- `\example` in EN-Version korrekt: `\newtheorem{example}` fehlt hier NICHT,
  da in paper26 keine `example`-Umgebung vorkommt ✓ (der Kommentar gilt nur für paper25)
- `\bibitem{TaylorWiles1995}` in EN-Version: wird durch `\cite` nicht referenziert —
  redundant, aber kein Fehler

**BUG-B6-07** (MINOR)
**Datei:** `paper26_l_function_elliptic_de.tex`, Zeile 230
**Problem:** Im historischen Beweis des Modularitätssatzes steht:
`Vermuted von Taniyama (1955) und Shimura (1957)`.
„Vermuted" ist kein deutsches Wort. Korrekt: „Vermutet" oder „Aufgestellt als Vermutung".

---

## Paper 27: Die Birch–Swinnerton-Dyer-Vermutung (EN + DE)

### Mathematische Korrektheit

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| Schwache BSD: ord_{s=1}L(E,s) = r | ✓ Korrekt | Als Vermutung korrekt formuliert |
| Starke BSD-Formel (eq:strong_bsd) | ✓ Korrekt | Alle Zutaten Ω_E, R_E, Sha, c_p, |tors|² korrekt |
| Tate-Shafarevich-Gruppe: alternierend ⟹ |Sha| Quadrat | ✓ Korrekt | |
| Coates-Wiles 1977 | ✓ Korrekt | CM-Bedingung korrekt |
| Kolyvagin 1990: r=0 bei L(E,1)≠0 | ✓ Korrekt | |
| Kolyvagin: r=1 bei L'(E,1)≠0 | ✓ Korrekt | |
| Gross-Zagier 1986: Heegner-Punkt ∞-Ordnung | ✓ Korrekt | |
| 2-Abstieg-Sequenz mit Selmer-Gruppe | ✓ Korrekt | |
| BSD-Produkt B_E(X) ~ C_E (log X)^r | ✓ Korrekt | Heuristische Begründung korrekt |

Die Starke BSD-Formel
$$\lim_{s\to 1} \frac{L(E,s)}{(s-1)^r} = \frac{\Omega_E \cdot |\Sha(E/\Q)| \cdot R_E \cdot \prod_p c_p}{|E(\Q)_{\mathrm{tors}}|^2}$$
ist korrekt — die übliche Formulierung mit reeller Periode Ω_E (nicht komplexe Periode).

### LaTeX-Korrektheit

**EN + DE:**
- `\Sha` als Operator korrekt deklariert ✓
- `\conjecture` in beiden Versionen definiert ✓
- Alle `\bibitem` vorhanden in beiden Versionen ✓

**BUG-B6-08** (MINOR)
**Datei:** `paper27_bsd_conjecture_de.tex`, Zeile 83
**Problem:** In der schwachen BSD-Vermutung steht `r = \mathrm{rang}\,E(\Q)`.
Im ganzen DE-Paper wird `\mathrm{rang}` verwendet (ohne Operator-Deklaration).
In der EN-Version ist es `\mathrm{rank}`. Der Operator `rang` sollte via
`\DeclareMathOperator{\rang}{rang}` deklariert oder durch `\operatorname{rang}`
ersetzt werden. Aktuell ist `\mathrm{rang}` syntaktisch korrekt, aber inkonsistent
mit dem Stil des Dokuments.

### EN↔DE Konsistenz

- Alle Inhalte konsistent ✓
- BSD-Formel identisch ✓

---

## Paper 28: Kongruente Zahlen, Elliptische Kurven und BSD (EN + DE)

### Mathematische Korrektheit

#### Theorem 2.1 (Äquivalenz): Kongruent ⟺ rank(E_n)≥1

**BUG-B6-09** (HOCH)
**Datei:** `paper28_congruent_numbers_bsd_en.tex`, Zeilen 117–165
**Problem:** Der Beweis der Richtung „Kongruent ⟹ Punkt unendlicher Ordnung"
ist unfertig und enthält sichtbare Entwurfsreste:
- Zeile 123: `y = ... (a suitable sign choice)` — Platzhalter!
- Zeile 128: `1/?` — Platzhalter!
- Zeile 151: `Hmm, let us use the direct formula instead.` — Entwurfskommentar!
- Zeile 157: `nein --` — Entwurfskommentar auf Deutsch im EN-Paper!
- Die Parametrisierung wird zweimal neu angefangen und der mittlere Versuch abgebrochen.

Der finale Ansatz (Zeilen 155–165) ist korrekt:
$(x_0,y_0) = ((c/2)², c(a²-b²)/8)$ liegt auf $E_n$, weil $c⁴-16n² = (a²-b²)²$.
ABER: Die Herleitung enthält den Fehler `c⁴-16n² = (a²-b²)²`. Prüfung:
- c²=a²+b², also c⁴=(a²+b²)²=a⁴+2a²b²+b⁴
- 16n²=16(ab/2)²=4a²b²
- c⁴-16n²=a⁴-2a²b²+b⁴=(a²-b²)² ✓

Der finale Beweis ist also mathematisch korrekt, aber die sichtbaren Entwurfsreste
(`Hmm`, `nein--`, `1/?`) machen das Paper als Publikation **nicht akzeptabel**.

**BUG-B6-10** (HOCH)
**Datei:** `paper28_congruent_numbers_bsd_en.tex`, Zeilen 287–302 (n=5, Beweis)
**Problem:** Der n=5-Beweis ist ebenfalls unfertig. Der Text beginnt die
Koordinatenberechnung, gibt dann `y_0 = 7/2 for the standard generator, but we can confirm...`
als unvollständige Aussage. Das Schluss-Statement lautet:
`we can confirm n=5 is congruent by the existence of the triangle (3/2, 20/3, 41/6)` —
was ausreichend wäre (da das Dreieck bereits im Example-Abschnitt steht), aber der
abgebrochene Berechnungsversuch davor ist ein Entwurfsrest.

**BUG-B6-11** (MITTEL)
**Datei:** `paper28_congruent_numbers_bsd_en.tex`, Zeilen 315–331
**Problem:** Numerischer Fehler in der Diskussion Tunnell-Kriterien für kleine n.
Der Text berechnet für n=3 (ungerade): A(3)=2, B(3)=1, also A(3)=2=2·1=2B(3).
Dann heißt es: „which would predict n=3 is congruent --- but this is the unresolved case".
Das ist korrekt: Unter BSD wäre n=3 kongruent, wenn A(3)=2B(3). ABER der Text sagt
zunächst (Zeile 319): „For n=2: C(2)=2, D(2)=1, so C(2)=2=2D(2), predicting 2 is congruent...
this contradicts expectation" — diese Selbstwiderlegung ist ein Entwurfsüberrest.
Die Werte A(1)=2, B(1)=1 stimmen NICHT: Für n=1 (ungerade) ist
2x²+y²+8z²=1 nur lösbar mit (x,y,z)=(0,±1,0): A(1)=2 ✓.
Für 2x²+y²+32z²=1: Nur (0,±1,0): B(1)=1 ✓.
Also A(1)=2=2·1=2B(1). Das Tunnell-Kriterium würde n=1 als kongruent vorhersagen.
ABER n=1 ist nicht kongruent (Fermat). Die Erklärung liegt im Fehler im Remark:
Das Kriterium A(n)=2B(n) ist notwendig für „n ist kongruent" (unconditional),
nicht hinreichend ohne BSD. Der Text erklärt dies korrekt weiter unten,
aber die Darstellung der Zählung ist verwirrend und inkonsistent.

**DE-Version (paper28_de.tex):**
Der entsprechende Beweisteil im DE-Paper ist besser, aber der gleiche erste Beweis-Fehler
(unklare Parametrisierung) wurde korrigiert: Die DE-Version startet direkt mit
`P = ((c/2)², c(a²-b²)/8)` ohne Entwurfskommentare.
Jedoch: Der abgebrochene Beweis auf der EN-Seite (Zeile 131-158: `nein --`, `Hmm`) ist
eine klare INKONSISTENZ zwischen EN und DE.

| Satz / Aussage | Status | Anmerkung |
|----------------|--------|-----------|
| n kongruent ⟺ rank(E_n)≥1 | ✓ Mathematisch korrekt | Beweis im EN-Paper hat Entwurfsreste |
| Torsion E_n(Q)_tors ≅ Z/2×Z/2 | ✓ Korrekt | Δ=64n⁶ ✓ |
| Formel (a,b,c) aus Punkt (x₀,y₀) | ✓ Korrekt | a²+b²=c² und ab/2=n verifiziert |
| Tunnell-Satz (Formulierung) | ✓ Korrekt | Zählformeln für A,B,C,D korrekt |
| n=6: Punkt (25/4, 35/8) | ✓ Korrekt | (35/8)²=1225/64=(25/4)³-36(25/4) ✓ |
| n=7: Dreieck (35/12, 24/5, 337/60) | ✓ Korrekt | Fläche und Hypotenuse verifiziert ✓ |

### LaTeX-Korrektheit

**EN-Version:**
- `\Fp` fehlt im Custom-Commands-Preamble — wird aber auch nicht direkt benutzt
- `\newtheorem{example}` korrekt definiert ✓
- Alle `\bibitem` vorhanden ✓

**DE-Version:**
- `\Fp` korrekt definiert ✓ — ABER wird in den Beweisen der DE-Version nicht verwendet ✓
- Alle `\bibitem` vorhanden ✓

---

## Zusammenfassung Batch 6

| Bug-ID | Schwere | Datei | Beschreibung |
|--------|---------|-------|--------------|
| BUG-B6-01 | MINOR | paper25_en.tex | j-Invariante doppelt/inkonsistent dargestellt |
| BUG-B6-02 | MINOR | paper25_en.tex | Singulärer Fall y²=x³-x²: b=0 Erklärung falsch |
| BUG-B6-03 | HOCH | paper25_en.tex | `example`-Umgebung nicht deklariert (LaTeX-Kompilierfehler) |
| BUG-B6-04 | MITTEL | paper25_de.tex | `\bibitem{Cornell1997}` fehlt in DE-Bibliographie |
| BUG-B6-05 | MITTEL | paper25_de.tex | Satz Néron-Ogg-Shafarevich fehlt komplett in DE-Version |
| BUG-B6-06 | MITTEL | paper26_en.tex | Entwurfsreste „(5,0)... wait," im Beispiel |
| BUG-B6-07 | MINOR | paper26_de.tex | Tippfehler: „Vermuted" statt „Vermutet" |
| BUG-B6-08 | MINOR | paper27_de.tex | `\mathrm{rang}` nicht als Operator deklariert |
| BUG-B6-09 | HOCH | paper28_en.tex | Beweis Thm 2.1: Entwurfsreste (`Hmm`, `nein--`, `1/?`) |
| BUG-B6-10 | HOCH | paper28_en.tex | n=5-Beweis unfertig: abgebrochene Berechnung |
| BUG-B6-11 | MITTEL | paper28_en.tex | Tunnell-Zahlen für kleine n: verwirrende/inkonsistente Darstellung |

### Prioritäten für Korrekturen

1. **Sofort beheben (HOCH):**
   - BUG-B6-03: `\newtheorem{example}[theorem]{Example}` in paper25_en.tex einfügen
   - BUG-B6-09 + BUG-B6-10: Beweisreste in paper28_en.tex entfernen und Beweis bereinigen

2. **Wichtig (MITTEL):**
   - BUG-B6-05: Néron-Ogg-Shafarevich in paper25_de.tex ergänzen
   - BUG-B6-06: Entwurfsrest in paper26_en.tex entfernen
   - BUG-B6-11: Tunnell-Diskussion in paper28_en.tex überarbeiten

3. **Kleinere Korrekturen (MINOR):**
   - BUG-B6-01, BUG-B6-02, BUG-B6-04, BUG-B6-07, BUG-B6-08

**Paper 28 (EN) muss vor Publikation überarbeitet werden.** Die mathematischen Inhalte
sind korrekt, aber der Beweis enthält sichtbare, nicht entfernte Entwurfskommentare,
die ein Publikationspaper disqualifizieren. Die DE-Version von Paper 28 ist deutlich
sauberer und kann als Vorlage für die EN-Bereinigung dienen.
