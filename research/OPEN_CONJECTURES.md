# OPEN_CONJECTURES.md - Offene mathematische Vermutungen
## Letzte Aktualisierung: 2026-03-08
## Autor: Kurt Ingwer

---

## Fernziel dieses Projekts
Das langfristige Ausbildungsziel ist das **Beweisen oder Widerlegen** offener
mathematischer Vermutungen. Dies erfordert den Aufbau tiefgreifenden Wissens
über Beweistechniken, Zahlentheorie, Analysis, algebraische Geometrie und mehr.

---

## I. Die Millennium-Preis-Probleme (Clay Mathematics Institute)

Das Clay Mathematics Institute hat 2000 sieben Probleme mit je 1 Million US-Dollar
Preisgeld ausgelobt. Stand 2026: Eines gelöst, sechs offen.

### 1. Riemann-Hypothese (RH) ⭐ PRIORITÄT 1
**Status:** Offen seit 1859
**Preisgeld:** 1.000.000 USD
**Formulierung:**
Alle nicht-trivialen Nullstellen der Riemann-Zeta-Funktion ζ(s) liegen auf der
"kritischen Geraden" Re(s) = 1/2.

Die Riemann-Zeta-Funktion:
```
ζ(s) = Σ(n=1 bis ∞) 1/n^s  (für Re(s) > 1)
```
durch analytische Fortsetzung auf ℂ \ {1} erweitert.

**Bekannte Fakten:**
- Über 10^13 Nullstellen numerisch auf der kritischen Geraden verifiziert
- Äquivalent zu Aussagen über die Verteilung der Primzahlen
- ζ(s) hat triviale Nullstellen bei s = -2, -4, -6, ...
- Funktionalgleichung: ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s)

**Verbindungen:**
- Primzahlverteilung (Primzahlsatz: π(x) ~ x/ln(x))
- L-Funktionen (GRH: Generalisierte Riemann-Hypothese)
- Zufallsmatrizen (GUE-Verteilung der Nullstellen-Abstände)
- Quantenchaos (Hilbert-Pólya-Vermutung: Nullstellen = Eigenwerte eines hermiteschen Operators)

**Ansätze für mich:**
- Numerische Verifikation weiterer Nullstellen
- Verbindung zu Quantenmechanik untersuchen
- Montgomery-Odlyzko-Gesetz der Nullstellen-Abstände
- Explizite Formel von Riemann: π(x) = Li(x) - Σ Li(x^ρ) + ...

---

### 2. P vs NP ⭐ PRIORITÄT 2
**Status:** Offen seit 1971 (Cook, Karp)
**Preisgeld:** 1.000.000 USD
**Formulierung:**
Ist jedes Problem, dessen Lösung in Polynomialzeit *verifiziert* werden kann,
auch in Polynomialzeit *lösbar*?

P = Klasse der in Polynomialzeit lösbaren Probleme
NP = Klasse der in Polynomialzeit verifizierbaren Probleme
Frage: P = NP?

**Bekannte Fakten:**
- Allgemeiner Konsens: P ≠ NP (aber kein Beweis)
- NP-vollständige Probleme: SAT, TSP, Clique, Vertex Cover, ...
- Wenn P = NP: Kryptographie bricht zusammen (RSA etc.)
- Baker-Gill-Solovay 1975: Orakel-Relativierungsmethoden reichen nicht aus
- Razborov-Rudich: "Natural Proofs" Barriere

**Ansätze für mich:**
- Komplexitätstheorie vertiefen (Schaltkreiskomplexität)
- Beweisbarrieren verstehen (Relativierung, Naturliche Beweise, Algebrisierung)
- Geometrische Komplexitätstheorie (Mulmuley-Sohoni)

---

### 3. Hodge-Vermutung
**Status:** Offen seit 1950
**Preisgeld:** 1.000.000 USD
**Formulierung:**
Auf einer nicht-singulären projektiven algebraischen Varietät über ℂ ist jede
Hodge-Klasse eine rationale Linearkombination von Klassen algebraischer Zyklen.

**Schwierigkeitsgrad:** Sehr hoch (erfordert algebraische Geometrie, Kohomologie)

---

### 4. Yang-Mills-Existenz und Masse-Gap
**Status:** Offen
**Preisgeld:** 1.000.000 USD
**Formulierung:**
Für jede kompakte einfache Eichgruppe G existiert eine nicht-triviale Quantisierung
der Yang-Mills-Theorie in ℝ⁴, und die zugehörige Hamilton-Operator-Theorie hat
eine Massenlücke Δ > 0.

**Verbindung:** Teilchenphysik (Quantenchromodynamik), starke Kernkraft

---

### 5. Navier-Stokes-Gleichungen ⭐ PRIORITÄT 3
**Status:** Offen
**Preisgeld:** 1.000.000 USD
**Formulierung:**
Existieren für die inkompressiblen Navier-Stokes-Gleichungen in ℝ³ glatte,
global definierte Lösungen für alle glatten Anfangsdaten mit endlicher Energie?

```
∂u/∂t + (u·∇)u = -∇p + νΔu + f
∇·u = 0
```

**Bekannte Fakten:**
- In 2D: Lösung bekannt (Ladyzhenskaya 1969)
- In 3D: Lokale Existenz bekannt, globale Regularität offen
- Verbindung zu Turbulenz

**Ansätze für mich:**
- Numerische Simulation (RK4, RK45 bereits implementiert)
- Reguläritätstheorie (Sobolev-Räume)
- Blow-up-Szenarien untersuchen

---

### 6. Birch-und-Swinnerton-Dyer-Vermutung
**Status:** Offen
**Preisgeld:** 1.000.000 USD
**Formulierung:**
Der Rang einer elliptischen Kurve E über ℚ ist gleich der Ordnung der Nullstelle
der L-Funktion L(E,s) bei s = 1.

**Verbindungen:** Kryptographie (elliptische Kurven-Kryptographie)

---

### 7. Poincaré-Vermutung ✅ GELÖST (2003)
**Status:** Bewiesen von Grigori Perelman
**Beweis-Methode:** Ricci-Fluss mit Chirurgie (Hamilton-Programm)
**Formulierung:** Jede einfach zusammenhängende, geschlossene 3-Mannigfaltigkeit
ist homöomorph zur 3-Sphäre.

**Lernwert:** Perelmans Methode (Ricci-Fluss) als Vorbild für geometrische Beweistechniken.

---

## II. Weitere bedeutende offene Vermutungen

### Goldbach-Vermutung ⭐ PRIORITÄT 1
**Status:** Offen seit 1742
**Formulierung:** Jede gerade ganze Zahl > 2 ist Summe zweier Primzahlen.
```
4=2+2, 6=3+3, 8=3+5, 10=3+7=5+5, ...
```
**Verifiziert bis:** 4 × 10^18
**Bekannt:** Ternäre Goldbach (Vinogradov 1937, vollständig bewiesen Helfgott 2013):
Jede ungerade Zahl > 5 ist Summe dreier Primzahlen.
**Ansatz:** Chen-Theorem (1966): Jede gerade Zahl = Primzahl + Produkt höchstens 2 Primzahlen

### Collatz-Vermutung (3n+1-Problem) ⭐ PRIORITÄT 1
**Status:** Offen seit 1937
**Formulierung:** Für jede positive ganze Zahl n führt die Folge
```
f(n) = n/2    (n gerade)
f(n) = 3n+1   (n ungerade)
```
nach endlich vielen Schritten zu 1.
**Verifiziert bis:** ~ 2^68 ≈ 2.95 × 10^20
**Terenz Tao (2019):** Fast alle Zahlen konvergieren (Dichte 1)
**Ansatz für mich:** Perfekt für numerische Exploration und Musteranalyse

### Zwillingsprimzahl-Vermutung
**Status:** Offen seit Antike
**Formulierung:** Es gibt unendlich viele Primzahlpaare (p, p+2).
```
(3,5), (5,7), (11,13), (17,19), (29,31), ...
```
**Bekannt:** Zhang 2013: Primpaare mit Abstand < 246 unendlich oft (→ Polymath-Projekt: 246)

### Legendre-Vermutung
**Formulierung:** Zwischen n² und (n+1)² liegt immer mindestens eine Primzahl.
**Folgt aus:** Cramér-Vermutung (gaps zwischen Primzahlen ≤ (ln p)²)

### ABC-Vermutung
**Formulierung:** Für coprime a,b,c mit a+b=c gilt fast immer c < rad(abc)^(1+ε)
**Status:** Shinichi Mochizuki behauptet Beweis (2012), IUT-Theorie noch umstritten

### Scholz-Vermutung
**Formulierung:** l(2^n - 1) ≤ n - 1 + l(n) (Additionsketten-Länge)

### Carmichael-Vermutung
**Formulierung:** Kein Wert der Eulerschen Phi-Funktion φ(n) tritt genau einmal auf.

---

## III. Strategischer Lernpfad zum Beweisen

### Phase 1: Werkzeuge (aktuell)
- [x] Algebra, Analysis, Lineare Algebra, Statistik, ODE
- [ ] Zahlentheorie (Primzahlsieb, Miller-Rabin, CRT)
- [ ] Komplexe Analysis (Riemann-Zeta, Residuensatz)
- [ ] Formale Beweistechniken (proof_theory.py)

### Phase 2: Grundlagen offener Probleme
- [ ] Analytische Zahlentheorie (für RH, Goldbach)
- [ ] Algebraische Geometrie (für Hodge, BSD)
- [ ] Komplexitätstheorie (für P vs NP)
- [ ] PDE-Theorie (für Navier-Stokes)

### Phase 3: Tiefenforschung
- [ ] Riemann-Zeta numerisch erkunden
- [ ] Goldbach/Collatz-Muster analysieren
- [ ] Verbindungen zwischen Problemen finden
- [ ] Teilresultate und Schranken beweisen

### Phase 4: Beweisversuche
- [ ] Collatz-Vermutung (als Einstieg)
- [ ] Goldbach für spezielle Klassen
- [ ] Beiträge zu Teilproblemen

---

## IV. Beweis-Techniken-Repertoire

| Technik | Anwendung | Beispiele |
|---------|-----------|-----------|
| Direkter Beweis | Implikation A→B | Pythagoras |
| Widerspruchsbeweis | Annahme führt zu Widerspruch | √2 irrational |
| Vollständige Induktion | Aussagen über ℕ | Gaußsche Summenformel |
| Starke Induktion | Komplexere ℕ-Aussagen | Fibonacci |
| Konstruktiver Beweis | Objekt explizit angeben | Primzahlexistenz |
| Diagonalargument | Überabzählbarkeit | Cantors ℝ > ℕ |
| Extremalprinzip | Kleinstes/größtes Element | Wohlordnung |
| Schubfachprinzip | Pigeonhole | Geburtstage |
| Probabilistischer Beweis | Zufälliges Objekt mit Eigenschaft | Ramsey-Theorie |
| Analytische Fortsetzung | Komplexe Funktionen | Riemann-Zeta |
| Algebraische Topologie | Invarianten berechnen | Poincaré |
| Modulformen | Zahlentheorie via Analysis | Fermat (Wiles) |

---

## V. Verbindungsnetz der Probleme

```
Primzahlen ──→ Riemann-Hypothese ──→ Goldbach (über Explizite Formel)
    │                │
    ▼                ▼
Miller-Rabin    L-Funktionen ──→ BSD-Vermutung
    │
    ▼
RSA-Kryptographie ──→ P vs NP (wenn P=NP: RSA bricht zusammen)

Navier-Stokes ──→ Turbulenz ──→ Yang-Mills (Physik-Verbindung)

Algebraische Geometrie ──→ Hodge ──→ Motiv-Theorie
                              │
                              ▼
                         BSD-Vermutung
```

---

## Quellen & Literatur

- Clay Mathematics Institute: https://www.claymath.org/millennium-problems/
- Terence Tao: "Structure and Randomness" (Aufsätze zu Primzahlen)
- Riemann (1859): "Über die Anzahl der Primzahlen unter einer gegebenen Größe"
- Perelman (2002-2003): Arxiv-Preprints zu Ricci-Fluss
- Andrew Wiles: Beweis des Großen Fermatschen Satzes (Modulformen-Technik)
- Timothy Gowers: "Mathematics: A Very Short Introduction"
