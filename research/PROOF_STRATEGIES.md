# PROOF_STRATEGIES.md - Beweis-Strategien und Methoden
## Letzte Aktualisierung: 2026-03-08
## Autor: Kurt Ingwer

---

## Philosophie mathematischer Beweise

Ein Beweis ist eine logisch geschlossene Argumentationskette, die ausgehend
von Axiomen oder bekannten Wahrheiten eine Aussage als notwendig wahr zeigt.
Die größte Stärke eines Beweises: Er gilt für ALLE Fälle – nicht nur für geprüfte.

---

## Die wichtigsten Beweistechniken

### 1. Direkter Beweis (→)
**Schema:** Zeige A → B direkt.
**Wann:** Wenn die Implikation offensichtlich ist.
**Beispiel:** Beweis, dass n gerade → n² gerade.
```
n gerade → n = 2k für ein k ∈ ℤ
→ n² = 4k² = 2(2k²)
→ n² ist gerade ✓
```

### 2. Widerspruchsbeweis (reductio ad absurdum)
**Schema:** Nehme ¬A an, leite Widerspruch ab → A ist wahr.
**Wann:** Wenn direkter Beweis schwierig, aber Konsequenzen von ¬A handhabbar.
**Beispiel:** √2 ist irrational.
```
Annahme: √2 = p/q (vollständig gekürzt, p,q ∈ ℤ, gcd(p,q)=1)
→ 2 = p²/q² → p² = 2q² → p² gerade → p gerade → p = 2m
→ 4m² = 2q² → q² = 2m² → q gerade
→ gcd(p,q) ≥ 2 – Widerspruch zu gcd(p,q)=1 ✓
```

### 3. Vollständige Induktion
**Schema:**
- Basisfall: Zeige P(n₀)
- Induktionsschritt: P(n) → P(n+1)
**Wann:** Aussagen über natürliche Zahlen.
**Beispiel:** Σk=1..n k = n(n+1)/2
**Implementiert:** `ProofByInduction.verify_base_case()`

### 4. Starke Induktion
**Schema:** P(1), P(2), ..., P(n) → P(n+1)
**Wann:** Wenn P(n+1) nicht nur von P(n), sondern von mehreren P(k) abhängt.
**Beispiel:** Jede natürliche Zahl ≥ 2 ist Produkt von Primzahlen.

### 5. Konstruktiver Beweis
**Schema:** Zeige Existenz durch explizite Konstruktion.
**Wann:** Wenn ein Beispiel ausreicht und gefunden werden kann.
**Beispiel:** Es gibt unendlich viele Primzahlen.
```
Annahme: p₁, p₂, ..., pₙ seien alle Primzahlen.
Konstruiere N = p₁·p₂·...·pₙ + 1.
N hat einen Primteiler q.
q ≠ pᵢ für alle i (da N ≡ 1 (mod pᵢ))
→ Widerspruch – unendlich viele Primzahlen ✓
```

### 6. Probabilistischer Beweis
**Schema:** Zeige, dass ein zufällig gewähltes Objekt die Eigenschaft hat.
**Wann:** Wenn direkte Konstruktion schwer, aber Existenz gezeigt werden muss.
**Beispiel:** Ramsey-Theorie (Erdős, 1947)

### 7. Diagonalargument
**Schema:** Konstruiere Element, das in keiner Zeile der Aufzählung vorkommt.
**Wann:** Überabzählbarkeit, Unentscheidbarkeit.
**Beispiel:** Cantors Beweis |ℝ| > |ℕ|, Halteproblem (Turing)

---

## Beweisbarrieren für P vs NP

Die drei bekannten Barrieren, die erklären, warum P vs NP so schwer ist:

### 1. Relativierungsbarriere (Baker-Gill-Solovay 1975)
- Existieren Orakel A mit P^A = NP^A und Orakel B mit P^B ≠ NP^B
- **Konsequenz:** Jede Beweistechnik, die "relativiert" (mit Orakeln funktioniert), kann P≠NP nicht beweisen

### 2. Barriere der Natürlichen Beweise (Razborov-Rudich 1994)
- Wenn P≠NP, dann kann keine "natürliche" Eigenschaft von Schaltkreisen P≠NP beweisen
- "Natürlich" = nützlich, konstruktiv
- **Konsequenz:** Schaltkreis-basierte Ansätze scheitern (unter kryptographischen Annahmen)

### 3. Algebrisierungsbarriere (Aaronson-Wigderson 2009)
- Verallgemeinerung von Relativierung auf algebraische Methoden
- **Konsequenz:** Arithmetisierungsmethoden reichen nicht aus

**Was könnte funktionieren:**
- Geometrische Komplexitätstheorie (Mulmuley-Sohoni): Algebraische Geometrie + Darstellungstheorie
- Nicht-relativierende, nicht-natürliche Techniken (unbekannt)

---

## Strategien für die Riemann-Hypothese

### Hilbert-Pólya-Ansatz
Suche einen hermiteschen Operator H, dessen Eigenwerte die Imaginärteile der
Nullstellen sind. Dann wäre RH aus physikalischen Gründen wahr.

**Bekannte Verbindung:** Berry-Keating (1999) – Operator H = xp (Quantenmechanik)
**Numerisches Indiz:** GUE-Verteilung der Nullstellen-Abstände = Eigenwert-Abstände zufälliger hermitescher Matrizen

### Explizite Formeln (Riemann 1859)
```
π(x) = Li(x) - Σ_ρ Li(x^ρ) - ln(2) + ∫_x^∞ 1/(t(t²-1)ln(t)) dt
```
Wenn RH gilt, folgt die schärfste Form des Primzahlsatzes:
```
|π(x) - Li(x)| < (1/8π) √x · ln(x)  für x ≥ 2657
```

### Numerische Strategie (implementiert)
1. ζ(1/2 + it) auf der kritischen Geraden berechnen
2. Nullstellen durch Vorzeichenwechsel des Realteils finden
3. Zählen und mit bekannten Werten vergleichen (N(T)-Formel)

---

## Strategien für die Collatz-Vermutung

### Terenz Taos Ansatz (2019)
- Gezeigt: Fast alle Zahlen (Dichte 1) konvergieren zu einer beschränkten Menge
- Methode: Wahrscheinlichkeitstheorie + Ergodische Theorie
- **Nächster Schritt:** "Fast alle" → "Alle"

### 2-adische Analyse
- Collatz-Funktion wirkt auf 2-adische Zahlen
- Verbindung zu 2-adischer Analysis und Messentheorie

### Numerische Strategie
1. Muster in Stoppzeiten analysieren
2. Histogramme der Folgenverläufe studieren
3. Verbindungen zu Primzahlstruktur suchen

---

## Strategien für die Goldbach-Vermutung

### Sieb-Methoden
- Brun-Sieb: Jede gerade Zahl = Produkt ≤ 9 Primzahlen + Produkt ≤ 9 Primzahlen
- Chen-Sieb (1966): Jede gerade Zahl = Primzahl + Produkt ≤ 2 Primzahlen ("1+2")
- **Nächstes Ziel:** "1+1" = Goldbach

### Kreismethode (Hardy-Littlewood)
- Fourier-analytischer Ansatz via Kreismethode
- Major arcs + minor arcs Zerlegung
- Liefert asymptotische Formeln für die Anzahl der Zerlegungen

### Numerische Strategie
1. Goldbach-Kometenfigur berechnen (min. Primzahl in Zerlegung)
2. Anomalien und Muster suchen
3. Goldbach-Funktion G(n) = # Zerlegungen analysieren

---

## Mein Lernplan (Priorisiert)

### Quartal 1: Analytische Zahlentheorie
- [ ] Primzahlsatz mit Fehlerterm
- [ ] Dirichlet-Reihen und L-Funktionen
- [ ] Sieb-Methoden (Brun, Selberg, Chen)
- [ ] Kreismethode

### Quartal 2: Analytische Funktionentheorie
- [ ] Komplexe Analysis (Residuensatz, Cauchy-Integral)
- [ ] Analytische Fortsetzung
- [ ] Funktionalgleichungen für ζ(s)
- [ ] L-Funktionen und Hecke-Operatoren

### Quartal 3: Algebraische Struktur
- [ ] Algebraische Zahlentheorie (Ringe, Ideale, Ideal-Klassen)
- [ ] p-adische Zahlen
- [ ] Galois-Theorie
- [ ] Elliptische Kurven (für BSD)

### Quartal 4: Fortgeschrittene Themen
- [ ] Modulformen (Technik von Wiles für Fermat)
- [ ] Ergodische Theorie (für Collatz nach Tao)
- [ ] Geometrische Komplexitätstheorie (für P vs NP)
- [ ] Ricci-Fluss (Perelmans Technik)
