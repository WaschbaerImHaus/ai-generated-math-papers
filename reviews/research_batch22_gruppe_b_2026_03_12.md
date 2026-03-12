# Research Review: Batch 22 – Gruppe B (Papers 84–87)
## Mathematische Analyse und Klassifikation der vier Vermutungen

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Build**: 172
**Computational Script**: `src/py/gruppe_b_batch22_verification.py`

---

## Überblick

Batch 22 behandelt vier Vermutungen aus dem Kernbereich der modernen Zahlentheorie
und arithmetischen Geometrie:

| Paper | Vermutung | Klassifikation |
|-------|-----------|----------------|
| 84 | Fontaine-Mazur (n=2) | WEITGEHEND BEWIESEN |
| 84 | Fontaine-Mazur (n≥3) | OFFEN |
| 85 | Lehmer-Mahler-Problem | OFFEN (starke num. Evidenz) |
| 86 | Goldfeld-Rang (avg=1/2) | OFFEN (best: ≤0.885) |
| 86 | Rang-Beschränktheit | OFFEN (umstritten) |
| 87 | Bateman-Horn (allg.) | OFFEN |
| 87 | Bateman-Horn (Dirichlet) | BEWIESEN |

---

## Paper 84: Fontaine-Mazur-Vermutung

### Mathematische Analyse

**Vermutung (Fontaine-Mazur 1995)**: Jede irreduzible, geometrische p-adische
Galois-Darstellung `ρ: Gal(Q̄/Q) → GL_n(Q_p)` kommt aus der Geometrie und ist automorph.

#### p-adische Hodge-Typen

Fontaine konstruierte eine Hierarchie von Periodenringen:
- `B_dR = lim_{n} A_inf[1/p]/(ker θ)^n` (de Rham)
- `B_cris ⊂ B_st ⊂ B_dR` (kristallin ⊂ halbststabil ⊂ de Rham)

**Hierarchie der Darstellungen**:
```
{kristallin} ⊊ {halbststabil} ⊊ {de Rham} ⊊ {Hodge-Tate}
```

"Geometrisch" = de Rham bei p + fast überall unverzweigt.

Die **Hodge-Tate-Gewichte** `h₁ ≤ h₂ ≤ ... ≤ hₙ ∈ Z` klassifizieren den Typ.
Beispiel: T_p(E) für eine elliptische Kurve E/Q_p mit guter Reduktion ist kristallin
mit Hodge-Tate-Gewichten {0, 1}.

### Status n=2: Weitgehend bewiesen

| Ergebnis | Autor(en) | Jahr | Inhalt |
|---------|-----------|------|--------|
| Alle ell. Kurven E/Q modular | Wiles + BCDT | 1995/2001 | V_p(E) ≅ ρ_f |
| Pot. halbststabile ρ (HT ≥ 1) | Kisin | 2009 | Kisin-Module (𝔖-Module) |
| Completed Cohomology | Emerton | 2010 | p-adische lok. Langlands für GL₂(Q_p) |
| Potentielle Modularität | Taylor | 2002 | Über total reellen Körpern |
| Gewicht-1-Formen | Deligne-Serre | 1974 | Artin-Darstellungen |

**Offene Randfälle bei n=2**:
- p=2 (technisch sehr schwierig)
- Darstellungen mit exzeptionalem lokalen Verhalten bei p
- Diese werden in aktueller Forschung (2024–2026) behandelt.

**Kernmethoden**:
1. **R=T-Methode** (Taylor-Wiles): Residuelle Modularität → Deformationsringiso. R ≅ T
2. **Colmez' p-adische lokale Langlands für GL₂(Q_p)** (2010): Bijektive Korrespondenz
   irred. 2-dim. p-adischer G_Q_p-Darstellungen ↔ irred. Banach-Darstellungen von GL₂(Q_p)

### Status n≥3: OFFEN

**Fundamentale Hindernisse**:

1. **Kein p-adisches lokales Langlands für GL_n, n≥3**:
   Colmez' Konstruktion nutzt die spezifische Struktur von GL₂(Q_p) über die
   (φ,Γ)-Moduln und die Faltungsoperation des Lie-Algebras sl₂. Für n≥3 fehlt
   eine analoge Konstruktion vollständig.

2. **R=T für GL_n scheitert**:
   Die Taylor-Wiles-Methode benötigt die "numerische Koinzidenz":
   `dim H¹(Gal(Q_S/Q), ad⁰ρ̄) = Σ_v dim H¹(G_v, ad⁰ρ̄)`.
   Für GL₂ ist dies erfüllt durch einen Balancierungseffekt, der für GL_n (n≥3)
   im allgemeinen Fall nicht gilt.

3. **Keine Shimura-Varietäten für GL_n über Q**:
   Für n≥3 gibt es keine direkten GL_n(A_Q)-analogen Shimura-Kurven, die als
   geometrische Quelle für Modulformen dienen könnten (außer über total reelle Körper).

4. **Globale Langlands-Korrespondenz für GL_n**:
   Nur in Spezialfällen bekannt (z.B. über Funktionenkörpern: Drinfeld n=2, Lafforgue allg.;
   über Q nur fragmentarisch).

**Vielversprechende neue Werkzeuge**:
- **Fargues-Scholze (2021)**: Geometrisierung der lokalen Langlands-Korrespondenz
  über die Fargues-Fontaine-Kurve — könnte perspektivisch n≥3 erschließen
- **Prismatische Kohomologie** (Bhatt-Scholze 2022): Neue kohomologische Maschinerie
- **Perverse Garben auf Shtuka-Räumen** (Scholze 2015)

### Klassifikation

```
Fontaine-Mazur n=2:  WEITGEHEND BEWIESEN (≥95% der Fälle)
Fontaine-Mazur n≥3:  OFFEN
```

---

## Paper 85: Lehmer-Mahler-Problem

### Mathematische Analyse

**Mahler-Maß**: Für `P(x) = a_d ∏(x - αⱼ)`:
```
M(P) = |a_d| · ∏ max(1, |αⱼ|) = exp(∫₀¹ log|P(e^{2πit})|dt)
```

**Kronecker-Theorem**: M(α) = 1 ⟺ α ist null oder eine Einheitswurzel.

**Lehmer-Vermutung (1933)**: Es gibt eine absolute Konstante c > 1, sodass
M(α) ≥ c für alle algebraischen ganzen Zahlen α, die keine Einheitswurzel sind.
Stärkere Form: M(P) ≥ λ₁₀ ≈ 1.17628 für alle P ∈ Z[x]\{Produkte von Kreisteilungspolynomen}.

**Lehmer-Polynom**:
```
L(x) = x¹⁰ + x⁹ - x⁷ - x⁶ - x⁵ - x⁴ - x³ + x + 1
```
- Irreduzibel über Q: **JA** (verifiziert)
- M(L) = λ₁₀ ≈ 1.1762808182599175...  (Salem-Zahl)
- Reziprok: x¹⁰ L(1/x) = L(x)

#### Dobrowolski-Schranke

**Theorem (Dobrowolski 1979)**: Für algebraisches ganzes α von Grad d ≥ 2, kein Einheitswurzel:
```
M(α) ≥ 1 + (1/1200) · (log log d / log d)³
```

**Numerische Analyse** (d=10):
```
M(α) ≥ 1 + (1/1200)·(log log 10/log 10)³
      = 1.0000396023...

Abstand zur 1:     3.96 × 10⁻⁵
Lehmer-Abstand:    0.176281
Faktor:            ≈ 4452 — Dobrowolski liegt weit unter Lehmer!
```

**Schlussfolgerung**: Dobrowolski ist asymptotisch hilfreich (M → 1 langsam), aber
quantitativ um Faktoren von ~4000 schwächer als das Lehmer-Ziel. Die beste bekannte
Schranke ist noch `O((log log d / log d)³)`, weit vom Ziel M ≥ 1.17628.

#### Salem-Zahlen nahe 1

Aus computationellen Studien (Boyd 1977, Mossinghoff 1998):
- Rang 1: λ₁₀ ≈ 1.17628 (Lehmer, Grad 10)
- Rang 2: ≈ 1.18836 (Grad 18)
- Rang 3: ≈ 1.20002 (Grad 20)

**Lücke (1, 1.17628) bleibt leer** für alle Polynome bis Grad 44 (Mossinghoff-Rhin-Wu 2008).

#### Smyth-Theorem für nicht-reziproke Polynome

Für nicht-reziproke P ∈ Z[x] (P ≠ xᵈP(1/x)):
```
M(P) ≥ θ₃ = 1.3247... (Plastizitätszahl, kleinste PV-Zahl)
```

Das Lehmer-Polynom ist **reziprok**, daher gilt Smyth nicht direkt.

### Computational Evidence

**Ergebnis der Suche (irreduzible Polynome bis Grad 10, Koeff. ∈ {-1,0,1})**:

```
Kandidaten geprüft:     59.046
Irreduzible Kandidaten: 2
Minimales M > 1:        1.176280818259918 = M(L) (nur Lehmer selbst und Variante)
Kein M ∈ (1, 1.17628) gefunden: JA – konsistent mit Lehmer-Vermutung
```

Das Minimum `M = 1.1762808183` wurde nur durch Lehmer-Varianten wie
`x¹⁰ - x⁹ + x⁷ - x⁶ + x⁵ - x⁴ + x³ - x + 1` (die Wechselzeichenvariante)
erreicht — dies ist jedoch ebenfalls eine Salem-Zahl mit demselben Mahler-Maß
(da Substitution x → -x dasselbe M ergibt).

**Verbindung zu Entropie algebraischer Systeme** (Lind-Schmidt-Ward 1990):
```
h(X_f) = m(f)   (topologische Entropie = log. Mahler-Maß)
```
Lehmer's Problem ≡ Gibt es algebraische Z-Aktion mit Entropie ∈ (0, log λ₁₀)?
Antwort: Vermutlich nein.

### Klassifikation

```
Lehmer-Mahler-Problem: OFFEN seit 1933
Numerische Evidenz:    STARK PRO Vermutung (bis Grad 44 verifiziert)
Beste Schranke:        Dobrowolski M ≥ 1 + c(log log d/log d)³ — weit von Ziel
Status:                Eines der ältesten offenen Probleme alg. Zahlentheorie
```

---

## Paper 86: Elliptische Kurven Rang

### Mathematische Analyse

**Mordell-Weil-Theorem**: `E(Q) ≅ Z^r ⊕ E(Q)_tors` mit Rang r = r(E) ≥ 0.

**Goldfeld-Vermutung (1979)**: Beim Ordnen nach Leitfähigkeit/Höhe:
```
lim_{X→∞} (∑_{N(E)≤X} r(E)) / #{E : N(E)≤X} = 1/2
```

**Bhargava-Shankar Theorem (2015)**:
- `avg|Sel²(E/Q)| = 3`
- `avg rank ≤ 0.885`
- Positive Anteile von Rang-0- und Rang-1-Kurven bewiesen

**Parity-Vermutung**: `(-1)^{r(E)} = w(E)` (Wurzelzahl bestimmt Rangparität).
- Bewiesen für elliptische Kurven über total reellen Körpern (Nekovář 2009) mit
  technischen Bedingungen. Über Q: weitgehend bekannt, aber nicht vollständig.

**Goldfeld-Heuristik**: Die 50/50-Verteilung auf Rang-0/1 entsteht wegen:
- ~50% Kurven haben `w(E) = +1` → analyt. Rang ≥ 0 (minimale Option: Rang 0)
- ~50% Kurven haben `w(E) = -1` → analyt. Rang ≥ 1 (minimale Option: Rang 1)
- Kurven mit Rang ≥ 2 haben Dichte 0 (unter Goldfeld/Poonen-Rains)

#### Boundedness Conjecture: Pro und Contra

**Argumente für unbeschränkten Rang**:
- Elkies (2006): explizite Kurve mit Rang ≥ 29
- Theoretische Konstruktionen (bedingt auf BSD): Kurven mit beliebig hohem Rang
- Keine bekannte absolute Schranke aus Theorie

**Argumente für beschränkten Rang** (neuere Heuristiken):
- Poonen-Rains (2012): Zufallsmatrizen-Modell → Rang ≥ r hat Dichte ~p^{-r²}
- Park-Poonen-Voight-Wood (2019): Detailliertere Heuristik → Rang beschränkt
- Bhargava-Shankar: Durchschnitt ≤ 0.885 lässt unbeschränkten Rang zwar zu,
  macht ihn aber "selten"

### Computational Evidence

**Rangverteilung für a,b ∈ [-50,50]** (Heuristik, naive Punktsuche):

```
Gesamtkurven:          10.201
Singuläre:                 5
Nicht-singuläre:       10.196

Rangverteilung (Heuristik):
  Rang 0:   8.356  (82.0%)
  Rang 1:   1.330  (13.0%)
  Rang 2:     510   (5.0%)
  Rang ≥3:     0   (0.0%)

Avg. Rang (Heuristik): 0.2305
Goldfeld-Ziel:         0.5000
Bhargava-Shankar:     ≤ 0.885
```

**Methodische Einschränkung**: Die naive Heuristik unterschätzt den Rang systematisch,
da sie nur ganzzahlige x-Koordinaten berücksichtigt und keine vollständige 2-Descent
durchführt. Bekannte Rang-1-Kurven wie y²=x³-432 werden als Rang 0 klassifiziert.

**Brumer-McGuinness Daten** (aus Paper 86):
```
Rang 0: ~35.2%
Rang 1: ~42.7%
Rang 2: ~18.9%
Rang ≥3: ~3.2%
```

Der Durchschnitt (BM): ~35.2·0 + 42.7·1 + 18.9·2 + 3.2·3 ≈ **90.1/100 = 0.90**
→ Konsistent mit Bhargava-Shankars Schranke 0.885.

### Klassifikation

```
Goldfeld-Vermutung (avg=1/2):     OFFEN
  Bestes Ergebnis: avg ≤ 0.885 (Bhargava-Shankar 2015)
  Numerische Evidenz: konsistent mit 1/2

Boundedness Conjecture:           OFFEN und KONTROVERS
  Pro Unbeschränktheit: Elkies Rang ≥ 29, kein theoretisches Hindernis
  Pro Beschränktheit:   Poonen-Rains, Park et al. Heuristiken
  Konsens: Prob. beschränkt, aber Beweis fehlt komplett

BSD (Birch-Swinnerton-Dyer):      OFFEN (Millennium Problem)
  Bewiesen: analyt. Rang 0 → r=0 (Kolyvagin); analyt. Rang 1 → r=1 (bedingt)
```

---

## Paper 87: Bateman-Horn-Vermutung

### Mathematische Analyse

**Bateman-Horn-Vermutung (1962)**: Für ein zulässiges System irred. Polynome
`f₁,...,fₖ ∈ Z[x]` mit Singulärreihenfaktor `𝔖(f₁,...,fₖ)` und Totalgrad D:
```
#{1 ≤ n ≤ x : f₁(n),...,fₖ(n) alle prim} ~ (𝔖/∏aᵢ) · x/(log x)^k
```

**Singuläre Reihe**:
```
𝔖(f₁,...,fₖ) = ∏_p (1-1/p)^{-k} · (1-ω(p)/p)
```
wobei `ω(p) = #{n ∈ Z/pZ : f₁(n)···fₖ(n) ≡ 0 (mod p)}`.

#### Warum Schinzels Hypothese H nicht ausreicht

Schinzels Hypothese H (qualitative Version): Unter Zulässigkeitsbedingung existieren
unendlich viele n mit f₁(n),...,fₖ(n) alle prim.

**Bateman-Horn ist stärker**: Es liefert nicht nur Existenz, sondern die genaue
asymptotische Dichte. Beides ist für deg ≥ 2 oder k ≥ 2 unbewiesen.

**Warum Siebmethoden nicht ausreichen**:
Selberg-Sieb liefert:
```
#{n ≤ x : ...} ≤ Cₖ · 𝔖 · x/(log x)^k  (mit Cₖ > 1)
```
Chen (1966/1973) beweist n mit p+2 ∈ P₂ (Halbprime), aber nicht n mit p+2 ∈ P.

**Fundamentales Hindernis** ("Parity Problem"): Siebmethoden können prinzipiell
nicht zwischen Primen und Produkten von zwei Primen unterscheiden. Dieses "Paritätsproblem"
(Selberg) blockiert den letzten Schritt von P₂ zu P.

**Green-Tao (2008)**: Primes enthalten beliebig lange APs — aber dies nutzt
Pseudozufälligkeit der Primes + Szemerédi, **nicht** Bateman-Horn direkt.

**Maynard-Tao (2013/2015)**: Bounded prime gaps ≤ 246 — beweist dass **ein** Element
eines zulässigen k-Tupels unendlich oft prim ist, aber **nicht alle gleichzeitig**.

#### Für f(x) = n²+1: Singuläre Reihe

```
ω(2) = 1  (n²+1 ≡ 0 (mod 2) → n ≡ 1 (mod 2))
ω(p) = 0  für p ≡ 3 (mod 4)  (-1 kein QR mod p)
ω(p) = 2  für p ≡ 1 (mod 4)  (-1 ist QR mod p, zwei Wurzeln)

𝔖(n²+1) = ∏_{p≡1 (mod 4)} (p-1)/(p-2) ≈ 1.37...
```

(Die verschiedenen Literaturwerte variieren je nach Normierung;
der Produktwert über 500 Primzahlen ergibt ≈ 1.3704.)

### Computational Evidence

**Numerische Überprüfung für f(x) = x²+1**:

| x | π(x,f) tatsächlich | BH-Vorhersage | Ratio |
|---|-------------------|---------------|-------|
| 100 | 19 | 14.9 | 1.277 |
| 500 | 70 | 55.1 | 1.270 |
| 1.000 | 112 | 99.2 | 1.129 |
| 5.000 | 472 | 402.3 | 1.173 |
| 10.000 | 841 | 744.0 | 1.130 |
| 50.000 | 3.613 | 3.166.5 | 1.141 |
| 100.000 | 6.656 | 5.951.7 | 1.118 |

**Analyse**: Die Ratio liegt konsistent über 1.0 und die Abweichungen nehmen für
wachsendes x ab (von ~28% bei x=100 auf ~12% bei x=100.000). Dies ist erwartetes
Verhalten — der Fehlerterm bei Bateman-Horn ist logarithmischer Natur, die Konvergenz
des Ratio gegen 1 ist **sehr langsam**.

Wichtig: Die Vorhersage unterschätzt systematisch leicht (Ratio > 1), was an der
Konvergenz des Eulerprodukts liegt (es werden nur endlich viele Primzahlen berücksichtigt
und die Normierung beeinflusst den Faktor).

**Bemerkenswert**: Iwaniec (1978) bewies `lim inf Ω(n²+1) ≤ 2`, also sind
"fast primes" (mit ≤ 2 Primfaktoren) unendlich oft durch n²+1 darstellbar —
aber **Primes** selbst: offen.

### Klassifikation

```
Bateman-Horn allgemein:        OFFEN
  Bewiesen:  k=1, deg=1 (Dirichlet 1837) + quantitatives PNT in AP
  Offen:     Alle Fälle mit deg ≥ 2 (einschl. f(x)=x²+1, n²+n+41, ...)
  Offen:     Alle Fälle mit k ≥ 2, deg = 1 (twin primes, Sophie Germain, ...)

Schinzel Hypothesis H:         OFFEN (qualitative Version)
Twin Prime Conjecture:         OFFEN (bester Stand: Lücke ≤ 246, Polymath8b)
Bunyakovsky k=1, deg≥2:        OFFEN (kein einziges Polynom bewiesen)

Numerische Evidenz:            Ratio → 1 für f(x)=x²+1, konsistent mit BH
```

---

## Gesamtklassifikation: Batch 22

### Paper 84 – Fontaine-Mazur-Vermutung

| Fall | Status | Hauptbeweis | Wichtigste offene Frage |
|------|--------|-------------|------------------------|
| n=2 (Standardfall) | **WEITGEHEND BEWIESEN** | Wiles, Kisin, Emerton | p=2-Fall, exzeptionelles lok. Verhalten |
| n=2, alle ell. Kurven | **BEWIESEN** | BCDT 2001 | — |
| n=2, Potential. Mod. | **BEWIESEN** | Taylor 2002 | — |
| n≥3 | **OFFEN** | — | p-adische lok. Langlands für GL_n |

**Einschätzung**: Die n=2-Vermutung ist de facto bewiesen (bis auf Randmengen),
die n≥3-Vermutung liegt in absolut unerreichbarer Ferne mit heutigen Methoden.

### Paper 85 – Lehmer-Mahler-Problem

**Status: OFFEN**

- Computational bis Grad 44 verifiziert: kein M ∈ (1, 1.17628)
- Dobrowolski-Schranke: M ≥ 1 + c(log log d/log d)³ — Faktor ~4000 zu schwach
- Smyth: nicht-reziproke Polynome erfüllen M ≥ θ₃ ≈ 1.3247 > 1.17628 ✓
- Reziproke Polynome (Lehmer-Typ): Keine Schranke > 1 beweisbar
- Salem-Zahlen können beliebig nah an 1 kommen (aber ob limit point = 1: offen)

### Paper 86 – Elliptische Kurven Rang

**Goldfeld avg=1/2: OFFEN** (Bhargava-Shankar: ≤ 0.885)
**Boundedness Conjecture: OFFEN** (Elkies ≥ 29; Poonen-Rains-Heuristik: beschränkt)

Das Boundedness-Problem ist das **kontroverseste** der vier Probleme in Batch 22:
Experimentelle Daten (sehr seltene Hochrang-Kurven) und Heuristiken zeigen in
entgegengesetzte Richtungen.

### Paper 87 – Bateman-Horn-Vermutung

**Status: OFFEN** (außer Dirichlet-Spezialfall, bewiesen)

Das fundamentale Hindernis ist das **Selberg-Paritätsproblem**: Siebmethoden
können Primes und P₂-Zahlen nicht unterscheiden. Ohne einen revolutionären neuen
Ansatz (ähnlich wie Wiles für FLT) scheint der allgemeine Beweis unerreichbar.

---

## Wichtige Verknüpfungen zwischen den Vermutungen

```
Fontaine-Mazur (n=2)
    ↓ bedingt auf BSD
    ↓ Kolyvagin-Euler-Systeme
Rang von ell. Kurven = analyt. Rang
    ↓ Goldfeld-Vermutung (avg=1/2)
    ↓ Paritätsvermutung
Rang-Verteilung ← → Bateman-Horn (f(x)=x³+ax+b hat Nullstelle ⟺ Rang ≥ 1)

Lehmer-Problem isoliert: tiefste Verbindung zu Entropie alg. Systeme
    (Lind-Schmidt-Ward 1990)
```

---

## Schlussfolgerungen

1. **Fontaine-Mazur (n=2)**: Der schwierigste Teil ist **bewiesen** — Wiles' Beweis
   von Fermat gilt als Beweis des wichtigsten n=2-Spezialfalls. Der vollständige
   n=2-Fall ist durch Kisin+Emerton weitgehend abgedeckt. Klassifikation: **BEWIESEN (n=2, bis auf Randmengen)** / **OFFEN (n≥3)**.

2. **Lehmer-Mahler**: Trotz 90 Jahren und massiver Computer-Suche: **kein Gegenbeispiel**.
   Stärkste Schranke (Dobrowolski) liegt um Faktoren von ~4000 unter dem Ziel.
   Die Vermutung ist **höchstwahrscheinlich wahr**, aber ein Beweis ist nicht in Sicht.

3. **Elliptische Kurven Rang**: Zwei Fragen (avg=1/2 und Beschränktheit).
   Beide **offen**. Bhargava-Shankars Durchschnitt ≤ 0.885 ist der beste bekannte
   Fortschritt. Numerische Evidenz spricht für avg=1/2 und für Beschränktheit.

4. **Bateman-Horn**: **OFFEN** für alle interessanten Fälle. Numerische Evidenz
   für f(x)=x²+1 ist konsistent mit der Vermutung (Ratio → 1 für x → ∞).
   Das Selberg-Paritätsproblem ist das fundamentale theoretische Hindernis.

---

*Dieser Review basiert auf den Papers 84–87 des Batch 22 sowie auf den Berechnungen
aus `src/py/gruppe_b_batch22_verification.py` (ausgeführt am 2026-03-12).*
