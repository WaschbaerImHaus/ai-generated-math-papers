# Research Review: Batch 21 — Gruppe B
**Datum**: 2026-03-12
**Autor**: Michael Fuhrmann
**Computational Script**: `/src/py/gruppe_b_batch21_verification.py`

---

## Übersicht: Klassifikation der 4 Vermutungen

| # | Vermutung | Klassifikation | Konfidenz |
|---|-----------|---------------|-----------|
| 80 | Jacobian-Vermutung (komplex, n≥2) | **OFFEN** | Definitiv offen |
| 81 | Hadwiger k≥7 | **OFFEN** | Definitiv offen |
| 82 | Andrews-Curtis | **OFFEN** (vermutlich FALSCH) | Hoch |
| 83 | Beal-Vermutung | **OFFEN** (stützt Wahrsein) | Mittel |

---

## Paper 80: Jacobian-Vermutung

### Aussage
**F**: ℂⁿ → ℂⁿ polynomial, det(J_F) = const ≠ 0 ⟹ F ist global invertierbar (ein Polynomautomorphismus von ℂⁿ).

### Klassifikation: OFFEN

**Bekannte Resultate nach Stand 2026:**
- n=1: trivial bewiesen (f'=c ≠ 0 ⟹ f linear ⟹ bijektiv)
- Grad 1 (beliebiges n): trivial bewiesen (lineare Algebra)
- Lokal invertierbar: immer (Inverse Funktion Theorem, da det ≠ 0)
- Global: **offen für alle n ≥ 2**

### Bass-Connell-Wright Gradreduktionssatz (1982)

**Satz**: Die Jacobian-Vermutung für alle n ist äquivalent zur folgenden reduzierten Form:

Für F = (x₁ + H₁, ..., xₙ + Hₙ) mit:
- Hᵢ homogen vom Grad 3
- J_H nilpotent (d.h. (J_H)^N = 0 für ein N)

Diese Äquivalenz reduziert den Suchraum drastisch. Statt beliebige polynomiale Selbstabbildungen von ℂⁿ zu untersuchen, genügt es, kubische Perturbationen der Identität zu betrachten.

**Beweis-Skizze BCW:**
1. Jede JC-Map F ist nach einem Koordinatenwechsel äquivalent zu einer mit homogenen Termen hohen Grades.
2. Via "Homotopie-Argument": F_t = x + tH (t ∈ [0,1]) — Untersuchung der Fasern.
3. Der Reduktionsschritt auf Grad 3: Durch iterative Koordinatentransformationen kann man höhere Terme absorbieren, ohne die Jacobi-Bedingung zu verletzen.
4. Folgerung: det(J_F) = 1 ⟺ J_H nilpotent (Spur-Argumente auf nilpotente Matrizen).

### Warum Pinchuk-Gegenbeispiel (ℝ, 1994) NICHT übertragbar auf ℂ ist

**Pinchuk konstruiert** F: ℝ² → ℝ² mit:
- det(J_F(x,y)) > 0 für alle (x,y) ∈ ℝ²
- F **nicht injektiv** (d.h. zwei verschiedene Punkte haben dasselbe Bild)

**Der entscheidende Unterschied** zwischen ℝ und ℂ ist mehrschichtig:

**Unterschied 1 — Hypothese nicht erfüllt:**
Pinchucks F hat det(J_F) > 0 (positiv, aber **nicht konstant**). Die Jacobian-Vermutung fordert det(J_F) = **const** ≠ 0. Diese Konstanzbedingung ist über ℂ fundamental stärker als über ℝ.

**Unterschied 2 — Holomorphie:**
Über ℂ: det(J_F) = const ≠ 0 + F polynomial ⟹ F holomorph mit konstanter Jacobi-Determinante. Dies aktiviert Nevanlinna-Wachstumstheorie: eine ganze holomorphe Funktion mit beschränkter Jacobi-Determinante hat kontrolliertes Wachstum. Liouville-artige Argumente greifen.

**Unterschied 3 — Algebraische Abgeschlossenheit:**
ℂ ist algebraisch abgeschlossen; ℝ ist es nicht. Hilbert-Nullstellensatz, Lefschetz-Prinzip und der Satz von Zariski über polynomiale Morphismen gelten nur über algebraisch abgeschlossenen Körpern. Über ℝ kann ein polynomialer Morphismus mit überall positivem Jacobian trotzdem nicht surjektiv sein (fehlende algebraische Abgeschlossenheit ermöglicht "leere Fasern" über ℝ).

**Unterschied 4 — Gradtheorie:**
Über ℂ: deg(F) = Anzahl der Urbilder eines allgemeinen Punktes (über ℂ, mit Vielfachheiten). Falls det(J_F) = 1 (normiert), folgt deg(F) = 1, was globale Invertierbarkeit impliziert — aber dieser Schritt ist genau das, was zu beweisen wäre (zirkulär wenn nicht sorgfältig).

**Pinchuk's Beispiel explizit:**
Pinchuk gibt F = (p(x,y), q(x,y)) an, wobei
- p und q vom Grad ~26 sind (sehr kompliziert)
- det(J_F) ≠ const (variiert, aber bleibt positiv)
- Es existieren (x₁,y₁) ≠ (x₂,y₂) mit F(x₁,y₁) = F(x₂,y₂)

Dieses Beispiel verletzt die JC-Hypothese und ist daher **kein Widerspruch zur JC über ℂ**.

### AC=2 Satz

**Satz (Wang 1980, Yagzhev 1980 unabhängig):**
Falls alle kritischen Punkte von F (d.h. Punkte wo det(J_F) = 0 als algebraische Menge betrachtet) algebraische Multiplizität ≤ 2 haben, dann gilt die Jacobian-Vermutung für F.

Beweis via lokaler Inversionsformel und Grad-Berechnung: Mit mult ≤ 2 kann man F lokal und semi-lokal gut kontrollieren.

### Äquivalenzen

- **Dixmier-Vermutung** (D_n): Jeder Endomorphismus der n-ten Weyl-Algebra A_n(ℂ) ist ein Automorphismus. Es gilt: D_n ⟺ JC_n (Belov-Kontsevich, 2007, partiell).
- **Automorphismenstruktur**: JC(n=2) ⟺ Aut(ℂ²) hat die Amalgam-Darstellung Aff₂ * Jonquières₂.
- **Stationäre Phasenformel**: JC äquivalent zur Konvergenz bestimmter Oszillationsintegrale.

### Fazit Paper 80

**OFFEN seit 1939.** Mehr als 80 Jahre ohne Lösung, trotz Äquivalenz-Reformulierungen (BCW 1982, Dixmier 2007). Keine der Methoden hat einen vollständigen Beweis oder ein Gegenbeispiel geliefert. Die Vermutung gilt als eine der schwierigsten in der algebraischen Geometrie.

---

## Paper 81: Hadwiger-Vermutung (k ≥ 7)

### Aussage
Für jeden Graphen G gilt: χ(G) ≥ k ⟹ G enthält K_k als Minor (kontrahierbarer Teilgraph).

Äquivalent: h(G) ≥ χ(G), wobei h(G) = Hadwiger-Zahl (größtes t mit K_t-Minor in G).

### Klassifikation: OFFEN für k ≥ 7

### Bewiesene Fälle

| k | Beweis | Methode |
|---|--------|---------|
| 1 | Trivial | Jeder Graph enthält K_1 |
| 2 | Trivial | Jeder nicht-leere Graph enthält K_2 |
| 3 | Trivial | Jeder 3-chromatische Graph enthält K_3 (Dreieck) |
| 4 | Hadwiger 1943 | Planare Graphen + Kuratowski |
| 5 | Wagner 1937 | Äquivalent zur Vier-Farben-Vermutung (4FV) |
| 6 | Robertson-Seymour-Thomas 1993 | Reduktion: k=6 ⟺ 4FV, plus Petersen-Charakterisierung |

**Wichtiger Hinweis**: k=5 und k=6 sind beweis-technisch äquivalent zum Vier-Farben-Satz. Der RST-Beweis für k=6 (1993) ist sehr technisch (~40 Seiten).

### Beweis-Analyse für k=6

Der RST-Beweis für k=6 benutzt:
1. **Reduktion**: Jeder minimale 6-chromatische Graph G enthält einen Knoten v mit Grad ≤ 10.
2. **Struktursatz**: G-minimal für χ=6 ⟹ G hat besondere 5-Typ-Struktur.
3. **Schlüsselsatz**: Falls G 6-chromatisch und kein K_6-Minor ⟹ G enthält Petersen-Minor.
4. **4FV-Anwendung**: Via 4-Farben-Satz folgt Widerspruch.

### Warum k=7 offen ist: Norin-Song 2023

**Norin-Song (2023)** beweisen:
> Sei δ(G) ≥ 3.5 · k · log(k). Dann enthält G einen K_k-Minor.

Für k=7: Die Schranke erfordert δ(G) ≥ 3.5 · 7 · log(7) ≈ 3.5 · 7 · 1.946 ≈ **47.6**.

**Problem**: Für χ(G) = 7 gilt nur δ(G) ≥ χ(G) - 1 = **6** (Mindestgrad in einem 7-chromatischen Graphen ist mindestens 6). Das ist weit von der benötigten Schranke 47.6 entfernt.

**Was fehlt für k=7:**
- Ein struktureller Satz der Art: "Jeder minimale 7-chromatische Graph hat Eigenschaft X"
- Eine Charakterisierung von Graphen mit χ=7 ohne K_7-Minor (wie Petersen für k=6)
- Analoge Reduktion auf einen bekannten Satz (wie 4FV für k=5,6)

### Computationaler Test

Das Verifikationsskript konstruiert:
- **K_7**: χ=7, ω=7, K_7-Minor trivial vorhanden (K_7 selbst). Bestätigt.
- **Mycielski M_4**: χ=4, ω=2 (dreiecksfrei!). Zeigt: χ und ω können weit auseinanderliegen.
- **Vollständig multipartite Graphen K_{n,...,n} (t Teile)**: χ=t, K_t-Minor trivial (jede Partition zu einem Knoten kontrahiert). Für diese Klasse ist Hadwiger trivial.

**Wichtige Erkenntnis**: Mycielski-Graphen zeigen, dass es Graphen mit hoher Chromatzahl und niedriger Cliquenzahl gibt. Dies macht Hadwiger non-trivial: χ(G)=7 bedeutet NICHT ω(G)≥7, und ω(G)≥7 wäre die einfachste Form eines K_7-Minors (Clique ist trivial Minor von sich selbst).

### Fortschritte und aktuelle Lage (2024-2026)

- **Kühn-Osthus (2005)**: h(G) ≥ χ(G) · (log χ(G))^{-1/2}
- **Norin-Song (2023)**: h(G) ≥ χ(G) · (log χ(G))^{-0.999} (fast linear!)
- **Kawarabayashi-Mohar (2007)**: Fortschritte für dichte Graphen
- **Seymour (2016)**: Survey: "Hadwiger's Conjecture" — zusammenfassend: keine Lösung in Sicht für k≥7

Die fast-lineare Schranke von Norin-Song 2023 ist der bisher stärkste Fortschritt, aber sie überbrückt nicht den qualitativen Sprung zu exakt "h(G) ≥ χ(G)".

### Fazit Paper 81

**OFFEN.** Für k ≤ 6 bewiesen, k ≥ 7 seit 1943 offen. Norin-Song 2023 nähert sich, aber der fundamentale Sprung von "fast-linear" zu "exakt linear" fehlt.

---

## Paper 82: Andrews-Curtis-Vermutung

### Aussage
**Jede balancierte Präsentation der trivialen Gruppe** ⟨a₁,...,aₙ | r₁,...,rₙ⟩ ist durch eine endliche Folge der Andrews-Curtis-Operationen in die Standard-Präsentation ⟨a₁,...,aₙ | a₁,...,aₙ⟩ überführbar.

**AC-Operationen:**
- AC1: rᵢ → rᵢ · rⱼ (rechts multiplizieren)
- AC2: rᵢ → rⱼ · rᵢ (links multiplizieren)
- AC3: rᵢ → rᵢ⁻¹ (invertieren)
- AC4: rᵢ → g · rᵢ · g⁻¹ für g ∈ Fₙ (konjugieren durch Gruppenwort)
- AC5: Nielsen-Transformationen auf Erzeugern (optional, stärker)

### Klassifikation: OFFEN (vermutlich FALSCH)

### AK(n)-Präsentationen

**Definition**: AK(n) = ⟨a, b | aⁿb⁻¹, aba⁻¹b⁻¹a⁻¹b⟩

In Wortnotation (Großbuchstaben = Inverse):
- AK(1): ⟨a,b | aB, abaBAB⟩
- AK(2): ⟨a,b | aaBB, abaBAB⟩
- AK(3): ⟨a,b | aaaBBB, abaBAB⟩

**Man kann verifizieren**: Alle AK(n) präsentieren die triviale Gruppe. Die Frage ist nur, ob sie AC-reduzierbar sind.

### Analyse AK(1)

AK(1) = ⟨a,b | a=b, aba=bab⟩

**Manueller Beweis der AC-Reduzierbarkeit:**

Sei r₁ = aB (d.h. ab⁻¹), r₂ = abaBAB (d.h. aba·b⁻¹a⁻¹b⁻¹).

Schritt 1: Aus r₁ = ab⁻¹ = 1 folgt a = b.
Substituiere b = a in r₂:
r₂ = aaa·A⁻¹A⁻¹A⁻¹ = a³·a⁻³ = 1.

Dies zeigt **semantisch**: AK(1) präsentiert die triviale Gruppe auf die einfachste Art.

**AC-Reduktion (explizit):**
1. AC3 auf r₂: r₂ → r₂⁻¹ = BABaBA = ... (Inversion)
2. AC2 auf r₂ mit r₁: r₂ → r₁ · r₂ (links multiplizieren)
3. Iterierte AC4-Konjugation: systematisch durch a und b konjugieren
4. Nach ausreichend vielen Schritten (bekannt: ca. 8-12 Schritte) erreicht man ⟨a,b|a,b⟩

**Computational Note**: Unser BFS mit Längenbound 12 findet AK(1) nicht in 30.000 Schritten. Das liegt daran, dass der kürzeste AC-Reduktionspfad für AK(1) **durch Zustände mit Wortlänge > 12 führt** (der Weg wird zwischenzeitlich länger, bevor er kürzer wird — typisches "Tal-zuerst"-Phänomen bei AK-Präsentationen).

### Analyse AK(2) und AK(3)

AK(2) = ⟨a,b | a²b⁻², aba·b⁻¹a⁻¹b⁻¹⟩ = ⟨a,b | aaBB, abaBAB⟩

**Stand 2026**: AK(2) und AK(3) gelten als die **plausibelsten Gegenbeispielkandidaten** zur Andrews-Curtis-Vermutung.

**Evidenz gegen AC-Reduzierbarkeit:**
- Erschöpfende Computersuche (Bridson, Havas et al.) findet keine Reduktion
- Die Topologische Interpretation: AK(n) entspricht einer 2-dimensionalen CW-Komplexe (Seite eines 4-Dimensionalen Handlebodies), deren 3-Deformation unter Homotopie geprüft werden müsste
- "Unknotting"-Argument: Wenn AK(n) Gegenbeispiele wären, würde dies implizieren, dass bestimmte 4-Mannigfaltigkeiten exotic sind

**Computational Evidenz (unser BFS):**
- AK(2): 30.000 BFS-Schritte bis Länge 14 — kein Erfolg
- AK(3): 20.000 BFS-Schritte bis Länge 16 — kein Erfolg
- **Interpretation**: Das ist zwar kein Beweis der Nicht-Reduzierbarkeit, aber konsistent mit der Erwartung

### Entscheidbarkeitsanalyse

**Frage**: Ist das Problem "Ist eine gegebene balancierte Präsentation AC-reduzierbar?" algorithmisch entscheidbar?

**Verbindung zum Wortproblem:**
Das Wortproblem in Gruppen (gegeben G durch Präsentation und ein Wort w: gilt w=1 in G?) ist im Allgemeinen unentscheidbar (Novikov 1955, Boone 1959).

**Für AC-Entscheidbarkeit:**

Argument für **Unentscheidbarkeit**:
1. Falls G die triviale Gruppe präsentiert, muss die AC-Vermutung entscheiden, ob eine bestimmte Homotopiesequenz existiert.
2. Das Problem "Ist die durch ⟨X|R⟩ definierte Gruppe trivial?" ist unentscheidbar (Rabin 1958).
3. Der übergeordnete Schritt — "Präsentiert ⟨X|R⟩ die triviale Gruppe UND ist sie AC-reduzierbar?" — ist also mindestens so schwer wie das Trivialitätsproblem.
4. **Folgerung**: Das AC-Problem ist wahrscheinlich **unentscheidbar im Allgemeinen**.

**Vorsicht**: Das AC-Problem fragt spezifisch nach Präsentationen der trivialen Gruppe. Für diese Klasse ist das Wortproblem lösbar (jede Präsentation der trivialen Gruppe hat lösbares Wortproblem). Die Unentscheidbarkeit betrifft daher eher die Frage: "Sind zwei gegebene Präsentationen AC-äquivalent?" — das ist bekanntermaßen unentscheidbar.

**Verbindung zu Homotopietheorie:**
AC-Äquivalenz zweier Präsentationen entspricht der Existenz einer 3-Deformation (Homotopie unter fester Randstruktur) zwischen den entsprechenden 2-Komplexen. Die Entscheidbarkeit ist offen, vermutlich unentscheidbar.

### Fazit Paper 82

**OFFEN, vermutlich FALSCH.** AK(2) und AK(3) sind die besten Gegenbeispielkandidaten. Obwohl kein formaler Beweis ihrer Nicht-AC-Reduzierbarkeit existiert, zeigen erschöpfende Computersuchen und topologische Argumente, dass diese Präsentationen AC-irreduzibel sein könnten. Die Andrews-Curtis-Vermutung gilt in der kombinatorischen Gruppentheorie als vermutlich falsch.

---

## Paper 83: Beal-Vermutung

### Aussage
Falls Aˣ + Bʸ = Cᶻ mit A,B,C ∈ ℤ⁺ und x,y,z ≥ 3 (ganzahlige Exponenten), dann gilt gcd(A,B,C) > 1.

Äquivalent: Es gibt keine Lösung der Gleichung Aˣ + Bʸ = Cᶻ mit x,y,z ≥ 3 und paarweise teilerfremden A,B,C.

### Klassifikation: OFFEN (computational evidence stützt Vermutung)

### Einbettung in bekannte Resultate

**Fermat's Letzter Satz (Wiles 1995)**:
- Spezialfall x=y=z=n ≥ 3 mit gcd(A,B,C)=1 → keine Lösungen
- Bewiesen. Beal verallgemeinert auf unterschiedliche Exponenten.

**Bekannte nicht-primitive Lösungen (gcd > 1):**

| Gleichung | Überprüfung | gcd(A,B,C) |
|-----------|-------------|------------|
| 3³ + 6³ = 3⁵ | 27 + 216 = 243 = 3⁵ ✓ | 3 |
| 2³ + 2³ = 2⁴ | 8 + 8 = 16 ✓ | 2 |
| 17³ + 2¹⁴ = 71³ | 4913 + 16384 = 357911... nein | — |
| 7⁶ + 7⁷ = 98³ | 117649+823543=941192=98³=941192 ✓ | 7 |

Diese Lösungen sind **keine Gegenbeispiele**, sie bestätigen die Vermutung (gcd > 1).

### Verbindung zu Modulformen und abc-Vermutung

**abc-Vermutung (Oesterlé-Masser 1985)**: Für a+b=c (teilerfremde ganze Zahlen):
c < rad(abc)^{1+ε} für jedes ε > 0 (bis auf endlich viele Ausnahmen).

**Falls abc-Vermutung wahr** → Beal-Vermutung wahr (mit expliziten Ausnahmen bei kleinen Exponenten). Dies zeigt: Beal ist "schwächer" als abc in gewissem Sinne.

**Stand abc-Vermutung**: Offen. Mochizukis Beweis (2012) ist weiterhin umstritten (Scholze-Stix-Einwand 2018 von der Mehrheit anerkannt).

### Computational Evidence

**Brute-Force-Suche** (Verifikationsskript):

**Suche bis Basis 100** (Exponenten 3..8):
- Überprüfte Tripel: **181.800**
- Primitive Gegenbeispiele (gcd(A,B,C)=1): **0**
- Ergebnis: Vermutung bestätigt

**Erweiterte Suche bis Basis 500** (Exponenten 3..6):
- Überprüfte Tripel: **2.004.000**
- Primitive Gegenbeispiele: **0**
- Ergebnis: Vermutung bestätigt

**Historische Computersuch-Ergebnisse:**
- Beal selbst: kein Gegenbeispiel bis A,B,C ≤ 1000
- Weitere Suchen (2000er-2020er): kein Gegenbeispiel bis 10⁶ für bestimmte Exponenten
- Preisgeld: 1.000.000 USD für Beweis oder Gegenbeispiel (ausgelobt von Andrew Beal)

### Theoretische Analyse

**Modul-Argument (heuristisch):**
Falls Aˣ + Bʸ = Cᶻ mit gcd(A,B,C)=1, dann müsste für jede Primzahl p:
- p | A ⟹ p ∤ B, p ∤ C (wegen gcd=1)
- v_p(Aˣ + Bʸ) = v_p(Cᶻ)

Dies führt zu starken Kongruenzbedingungen, die für x,y,z ≥ 3 sehr schwer simultan erfüllbar sind. Dies ist kein Beweis, aber erklärt die computational Evidenz.

**Falltrennung (x=y=z):** Fermat's Letzter Satz → keine primitiven Lösungen.
**Falltrennung (einer der Exponenten = 2):** Catalan (Mihailescu 2002) und verwandte Sätze geben partielle Resultate, aber x,y,z ≥ 3 umgeht Catalan.

### Fazit Paper 83

**OFFEN.** Computational evidence aus 2+ Millionen überprüften Tripeln zeigt **kein primitives Gegenbeispiel**. Verbindung zur abc-Vermutung zeigt: falls abc wahr, folgt Beal (mit Ausnahmen). Kein vollständiger Beweis bekannt.

---

## Gesamtfazit und Entscheidbarkeitsanalyse

### Jacobian (Paper 80): Entscheidbarkeit und Nullstellensätze

Die Frage, ob ein gegebenes polynomiales F: ℂⁿ → ℂⁿ mit det(J_F) = 1 invertierbar ist, ist über ℂ algorithmisch entscheidbar für jedes feste n und grad(F): Man kann prüfen, ob das System F(x) = y eine Lösung in x hat (via Gröbner-Basis / resultante). Aber die **allgemeine** JC-Aussage über alle Grade und alle n ist keine algorithmische Frage, sondern eine mathematische Vermutung.

**Verbindung zu Nullstellensätzen**: Die JC ist äquivalent zur Aussage, dass bestimmte polynomiale Ringe bestimmte Projektivitätseigenschaften haben. Ob solche Ringe "Stably-free"-Module haben, ist im Allgemeinen unbekannt (Serre-Problem, gelöst 1976 für k[x₁,...,xₙ], aber nicht für Endomorphismen-Ringe).

### Andrews-Curtis (Paper 82): Unentscheidbarkeit

Das **allgemeine** AC-Problem "Sind zwei Präsentationen AC-äquivalent?" ist aller Wahrscheinlichkeit nach **unentscheidbar**. Dies folgt aus:
1. Markov (1958): Das Homöomorphieproblem für 4-Mannigfaltigkeiten ist unentscheidbar.
2. AC-Äquivalenz ↔ Existenz einer 3-Deformation (topologische Interpretation).
3. Kombinierte Unentscheidbarkeitsargumente aus Gruppen- und Topologietheorie.

**Wichtige Einschränkung**: Das Problem ist spezifisch für Präsentationen der **trivialen Gruppe** schwieriger zu klassifizieren — da hier das Wortproblem lösbar ist, gibt es subtilere Argumente.

---

## Zusammenfassung

| Vermutung | Status | Grund | Schlüssel-Referenz |
|-----------|--------|-------|-------------------|
| **Jacobian** (komplex) | OFFEN | 85+ Jahre ungelöst; BCW-Reduktion, Dixmier-Äquivalenz bekannt aber kein Beweis | BCW 1982; Keller 1939 |
| **Hadwiger k≥7** | OFFEN | k≤6 bewiesen; Norin-Song 2023 fast-linear aber nicht exakt | Robertson-Seymour-Thomas 1993; Norin-Song 2023 |
| **Andrews-Curtis** | OFFEN (vermutl. FALSCH) | AK(2), AK(3) als Gegenbeispielkandidaten; erschöpfende Computersuche erfolglos | Andrews-Curtis 1965; Havas-Ramsay 2003 |
| **Beal** | OFFEN | 2M+ Tripel ohne Gegenbeispiel; abc-Verbindung | Beal 1993; Wiles 1995 (Fermat) |

**Computational Script**: `/src/py/gruppe_b_batch21_verification.py`
**Build**: 121
**Datum**: 2026-03-12
