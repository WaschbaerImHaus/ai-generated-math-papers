# Vollstaendiges Mathematisches Audit -- Batch 1 & Batch 2 (Papers 1--12)

**Build 14 -- Datum: 2026-03-12**
**Reviewer: Claude Opus 4.6 (Peer-Review-Modus)**

---

## Zusammenfassung

24 Dateien geprueft (12 Papers, je EN + DE). Alle Beweise wurden Schritt fuer Schritt verifiziert, Formeln nachgerechnet, Definitionen auf Konsistenz geprueft, EN- und DE-Versionen verglichen.

**Gesamtbewertung:**
- 9 von 12 Papers: DRUCKREIF (keine oder nur geringfuegige Maengel)
- 2 Papers: UEBERARBEITUNG NOETIG (Paper 2 EN + DE)
- 1 Paper: DRUCKREIF mit geringer Anmerkung (Paper 8)

---

## MAENGELLISTE

### =====================================================================
### PAPER 1 -- Giuga Three Primes (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

Alle Beweise sind korrekt und lueckenlos:

1. **Lemma Squarefreeness (Zeilen 128-133 EN / 129-134 DE):** Korrekt. Wenn p^2 | n, dann n/p = 0 mod p, also n/p - 1 = -1 mod p, und p | -1 ist unmoeglich.

2. **Lemma Two-prime (Zeilen 140-148 EN):** Korrekt. q | (p-1) mit p < q erzwingt p-1 = 0.

3. **Theorem 2qr (Zeilen 158-187 EN):** Korrekt. r = 2q - 1 aus der schwachen Bedingung, dann Paritaetswiderspruch bei der starken Bedingung: 2(q-1) | (2q-1) ist unmoeglich da gerade | ungerade.

4. **Theorem pqr all odd (Zeilen 200-253 EN):** Korrekt.
   - r(r-1) | (pq-1) folgt korrekt aus gcd(r, r-1) = 1
   - r >= q+2 fuer aufeinanderfolgende ungerade Primzahlen: korrekt
   - r(r-1) >= (q+2)(q+1) = q^2 + 3q + 2: korrekt
   - pq >= q^2 + 3q + 3, also p >= q + 3 + 3/q: korrekt (nachgerechnet)
   - Da q >= 5 gilt 3/q < 1, also p >= q + 4: korrekt
   - Widerspruch zu p < q: korrekt

5. **Literatur (Zeile 87 EN):** "Borwein, Borwein, Girgensohn, and Parnes" -- KORREKT (fruehere Versionen hatten hier einen Fehler, jetzt behoben).

6. **Numerische Verifikation:** Giuga-Zahlen {30, 858, 1722, 66198} alle korrekt faktorisiert.

7. **EN vs DE Vergleich:** Inhaltlich identisch. Keine mathematisch relevanten Uebersetzungsabweichungen. DE verwendet (S)/(St) statt (W)/(S) als Labels -- konsistent innerhalb jeder Version.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine neuen Fehler gefunden |

---

### =====================================================================
### PAPER 2 -- Lehmer Three Primes (EN + DE)
### =====================================================================

**Urteil: UEBERARBEITUNG NOETIG**

#### Bekannte Bugs -- Verifikation:

| Bug-ID | Status | Verifikation |
|--------|--------|-------------|
| BUG-B1-P2-EN-001 | **BESTAETIGT** | Zeilen 310-314 EN: "q >= 17: One verifies that phi(5qr) = 4(q-1)(r-1) grows faster than 5qr - 1... An exhaustive computer-verified check for all prime q <= 500 confirms no valid triple (5, q, r) exists." -- Dies ist KEIN elementarer Beweis. Fuer p=5, q>=17 wird nur auf eine Computersuche verwiesen. Ein mathematisches Paper, das "elementary" im Titel fuehrt, sollte hier einen analytischen Beweis liefern. |
| BUG-B1-P2-EN-002 | **BESTAETIGT** | Zeilen 316-331 EN: Im Fall p>=7 wird die Schranke A >= (p-1)(q-1)/(p+q-2) > q*(p-2)/(p-2+1) >= q/2 benutzt. Die Ungleichung ist korrekt, aber der Uebergang von "A <= 2p+q-2" zu "q < 6(p-1)/(p-5)" ist unzureichend kommentiert und wird dem Leser nicht nachvollziehbar gemacht. |
| BUG-B1-P2-DE-001 | **BESTAETIGT** | Zeilen 216-219 DE: In Lemma schranke_k wird berechnet: k = (pq-1)/(r-1) <= (pq-1)/(q+1) = p - 1 - p/(q+1) + (p-1)/(q+1). Diese Zerlegung ist algebraisch FALSCH. Korrekt: (pq-1)/(q+1) = p - (p+1)/(q+1). Die Schlussfolgerung k <= p-1 ist trotzdem richtig (da (pq-1)/(q+1) < p), aber der Zwischenschritt ist fehlerhaft. |

#### Neue Befunde:

| ID | Paper | Sprache | Schwere | Zeile | Beschreibung | Fix |
|----|-------|---------|---------|-------|-------------|-----|
| BUG-B14-P2-EN-003 | 2 | EN | GERING | 132 | "since q >= p + 2 >= 5 we have q - 2 >= 3 >= 2" -- Die Aussage q >= p+2 setzt voraus, dass p und q verschiedene ungerade Primzahlen sind. Fuer p=3 folgt q >= 5, korrekt. Aber die Formulierung "q >= p + 2" ist nur fuer ungerade p+q korrekt; falls p=2, q=3 waere q = p+1. Der Text schraenkt aber korrekt auf p >= 3 ein (Case p >= 3). | Kein Fix noetig, nur Klarstellung wuenschenswert |
| BUG-B14-P2-EN-004 | 2 | EN | MITTEL | 248-252 | Die Reduktion von (p-1)(q-1) | (r-1)(alpha+k) zu "pq equiv alpha pmod{(p-1)(q-1)}" beruft sich auf eine Kongruenz, die nicht explizit bewiesen wird. Der Schritt pq = (p-1)(q-1) + (p-1) + (q-1) = (p-1)(q-1) + p+q-2 ist korrekt, aber der Uebergang pq+k equiv alpha+k benoetig einen klareren Beweis. | Zwischenschritt explizit auffuehren |
| BUG-B14-P2-DE-002 | 2 | DE | GERING | 250 | "(p-1) | (q-1) eine notwendige Vereinfachung liefert (aus n equiv 1 pmod{p-1} folgt qr equiv 1 pmod{p-1})" -- Dies ist NICHT korrekt als "notwendige Vereinfachung". Aus qr equiv 1 mod (p-1) folgt NICHT (p-1) | (q-1). Es folgt nur q equiv r^{-1} mod (p-1). Die Aussage (p-1) | (q-1) wird zwar spaeter fuer die Fallanalyse benutzt, ist aber hier falsch begruendet. | Begruendung korrigieren: (p-1) | (q-1) folgt nicht direkt aus qr equiv 1 mod (p-1) |
| BUG-B14-P2-DE-003 | 2 | DE | MITTEL | 316-318 | "Fuer p = 11: q < 6*10/6 = 10 < 11 = p" -- Die Rechnung 6*10/6 = 10 ist korrekt, aber die Formel q < 6(p-1)/(p-5) wird hier mit p=11 verwendet: 6*10/(11-5) = 60/6 = 10. Da q > p = 11 verlangt wird, ist q > 11 > 10 ein Widerspruch. Korrekt. Aber "Fuer p = 13: q < 6*12/8 = 9 < 13" -- 6*12/8 = 9 ist korrekt, und q > 13 > 9. Korrekt. Aber die Formel wurde in Lemma grenze_q nur fuer D=1 bewiesen, nicht allgemein. Fuer D > 1 wird ohne Beweis behauptet, dass dieselbe Schranke gilt. | Lemma fuer allgemeines D beweisen oder D=1-Einschraenkung dokumentieren |

---

### =====================================================================
### PAPER 3 -- Giuga-Carmichael (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

1. **Theorem Congruence (Zeilen 126-148 EN):** Korrekt.
   - Aus Giuga (W): n/p equiv 1 mod p, also n equiv p mod p^2. Korrekt.
   - Aus Korselt: (p-1) | (n-1), also n equiv 1 mod (p-1), also n equiv p mod (p-1). Korrekt.
   - gcd(p^2, p-1) = 1 fuer Primzahlen p: Korrekt (jeder gemeinsame Teiler d von p^2 und p-1 teilt gcd(p, p-1)=1).
   - CRT gibt n equiv p mod p^2(p-1). Korrekt.

2. **CRT Obstruction (Zeilen 160-197 EN):** Nachgerechnet.
   - n equiv 3 mod 18, n equiv 5 mod 100, n equiv 7 mod 294.
   - n = 18k + 3, einsetzen: 18k + 3 equiv 5 mod 100, also 18k equiv 2 mod 100, also 9k equiv 1 mod 50.
   - 9 * 39 = 351 = 7*50 + 1. Korrekt: 9^{-1} equiv 39 mod 50.
   - k equiv 39 mod 50, also n = 900m + 705. Korrekt: lcm(18, 100) = 900.
   - gcd(900, 294): 900 = 3*294 + 18, 294 = 16*18 + 6, 18 = 3*6. Also gcd = 6. Korrekt.
   - 7 - 705 = -698. -698 mod 6: 698 = 116*6 + 2, also -698 equiv -2 equiv 4 mod 6. Korrekt. 4 != 0, also unlösbar.

3. **Proposition odd (Zeilen 106-116 EN):** Korrekt. (-1)^{n-1} = -1 fuer gerades n, also -1 equiv 1 mod n, also n | 2, Widerspruch zu n > 2 zusammengesetzt.

4. **Moduli korrekt:** p=3: 9*2=18. p=5: 25*4=100. p=7: 49*6=294. Alle korrekt.

5. **EN vs DE:** Inhaltlich identisch. Keine Abweichungen.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine Fehler gefunden |

---

### =====================================================================
### PAPER 4 -- Giuga Squarefree (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

Kurzes Paper mit einem einzigen zentralen Beweis. Mathematisch trivial und korrekt: p^2 | n => p | n/p => n/p equiv 0 mod p => n/p - 1 equiv -1 mod p => p | -1, Widerspruch.

Giuga-Zahlen {30, 858, 1722, 66198} korrekt faktorisiert. 66198 = 2 * 3 * 11 * 17 * 59. Nachgerechnet: 2*3*11*17*59 = 6*11*17*59 = 66*17*59 = 1122*59 = 66198. Korrekt.

EN vs DE: Identisch.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine Fehler gefunden |

---

### =====================================================================
### PAPER 5 -- Giuga No Semiprime (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

Beweis korrekt: Fuer n = pq mit p < q: schwache Bedingung (W) fuer q ergibt q | (p-1). Da 0 <= p-1 < q, muss p-1 = 0, also p = 1, Widerspruch.

EN vs DE: Identisch.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine Fehler gefunden |

---

### =====================================================================
### PAPER 6 -- Lehmer Squarefree (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

Beweis korrekt: p^2 | n => p | phi(n) (da p^{a-1}(p-1) | phi(n)) => p | (n-1) (da phi(n) | (n-1)) => p | n - (n-1) = 1. Widerspruch.

Verallgemeinerung auf M*phi(n) = n-1 korrekt angemerkt.

EN vs DE: Identisch.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine Fehler gefunden |

---

### =====================================================================
### PAPER 7 -- Lehmer No Semiprime (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF (DE mit Stilanmerkung)**

Beweis korrekt:
- Identitaet pq - 1 = (p-1)(q-1) + (p-1) + (q-1): Nachgerechnet, korrekt.
- Daher pq - 1 equiv p + q - 2 mod (p-1)(q-1).
- Fall p=2: (q-1) | q. Da q = (q-1)+1, folgt (q-1) | 1, also q = 2 = p. Widerspruch zu q > p.
- Fall p>=3: (p-2)(q-2) - 1 >= 1*2 - 1 = 1 > 0. Korrekt, da p >= 3 => p-2 >= 1 und q >= p+2 >= 5 => q-2 >= 3.

| ID | Paper | Sprache | Schwere | Zeile | Beschreibung | Fix |
|----|-------|---------|---------|-------|-------------|-----|
| BUG-011 | 7 | DE | GERING | -- | Stilanmerkung: Der DE-Text in Remark (Zeile 120-122) sagt "(p-2)(q-2) > 1, was fuer alle Primpaare p >= 3, q > p gilt (da p-2 >= 1 und q-2 >= p-1 >= 2)". Die Behauptung q-2 >= p-1 ist nur korrekt wenn q >= p+1, d.h. q >= p+2 (aufeinanderfolgende ungerade Primzahlen). Fuer den Fall p=2 ist q-2 = q-2 und p-1 = 1, also q-2 >= 1 ist trivial. Der Beweis behandelt p=2 separat, daher kein mathematischer Fehler, aber die Anmerkung ist leicht irrefuehrend. | Stilistische Ueberarbeitung |

---

### =====================================================================
### PAPER 8 -- Wilson's Theorem (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF (mit geringer Anmerkung)**

1. **Wilson if-direction (Zeilen 112-138 EN):** Korrekt. Paarungsargument: Selbst-inverse Elemente sind 1 und p-1 (Lemma korrekt). (p-3)/2 Paare, jedes Paar gibt 1. Produkt = 1 * 1^{(p-3)/2} * (p-1) = -1 mod p.

2. **Wilson only-if (Zeilen 142-160 EN):** Korrekt.
   - Fall a < b: n = ab | (n-1)!, da a,b verschiedene Faktoren in {1,...,n-1}.
   - Fall a = b, a >= 3: a und 2a in {1,...,n-1} (da 2a < a^2 fuer a >= 3), also 2n | (n-1)!.
   - Fall a = 2, n = 4: (4-1)! = 6 equiv 2 mod 4. Korrekt.

3. **Composite factorial (Zeilen 166-185 EN):** Fuer n >= 6 zusammengesetzt: (n-1)! equiv 0 mod n. Korrekt. Ausschluss von n=4 (a=2) explizit.

4. **Wilson for p^k (Zeilen 221-256 EN):** Korrekt.
   - Ungerades p: |S| = 2 (nur +/-1), also Produkt = -1.
   - p=2, k=1: Produkt = 1 equiv -1 mod 2. Korrekt.
   - p=2, k=2: {1,3}, Produkt = 3 equiv -1 mod 4. Korrekt.
   - p=2, k>=3: Vier Loesungen, Produkt = -(2^{2k-2}-1). Fuer k >= 3: 2k-2 >= 4 >= k, also 2^{2k-2} equiv 0 mod 2^k, also Produkt equiv 1 mod 2^k. Korrekt.

5. **General Wilson (Zeilen 300-335 EN):** Korrekt. Drei Faelle: |S|=1 => e, |S|=2 => tau, |S|>=4 => e.

| ID | Paper | Sprache | Schwere | Zeile | Beschreibung | Fix |
|----|-------|---------|---------|-------|-------------|-----|
| BUG-B2-P8-EN-001 | 8 | EN+DE | GERING | -- | Die Ungleichung "2a < a^2 fuer a >= 3" in der only-if-Richtung wird auch als "2a = 2sqrt(n) < n fuer n >= 9" geschrieben (Zeile 155 EN). Fuer k = a = 3 (n=9) gilt 2*3 = 6 < 9, korrekt. Aber die Behauptung "2k-2 >= k+1" fuer die 2^k-Analyse (Zeile 253 EN: "2k-2 >= k+1 > k") gilt NICHT fuer k=3: 2*3-2 = 4, k+1 = 4, also 4 >= 4 ist Gleichheit, nicht strikt. Die Aussage "2k-2 >= k" (die genuegt) gilt ab k >= 2, was ausreicht. Die unnoetig scharfe Behauptung "2k-2 >= k+1" ist fuer k=3 falsch (Gleichheit, nicht >). | Korrektur: "2k-2 >= k" statt "2k-2 >= k+1 > k" |

---

### =====================================================================
### PAPER 9 -- Wilson Prime Powers (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

1. **Lemma sqrt_unity (Zeilen 85-105 EN):** Korrekt. min(alpha, beta) = 0 da p ungerade und Differenz = 2.

2. **Theorem wilson_pk_odd (Zeilen 118-154 EN):** Korrekt. Paarung, |G|-2 ist gerade (p-1 gerade), Produkt = 1 * (-1) = -1.

3. **Theorem wilson_pk_2 (Zeilen 163-209 EN):** Korrekt. k=1: trivial. k=2: {1,3}, Produkt=3=-1. k>=3: 4 selbst-inverse Elemente, Produkt = -(2^{2k-2}-1) equiv 1 mod 2^k.

4. **Structure Theorem (Zeilen 218-260 EN):** Korrekt. Lifting-Argument fuer Primitivwurzeln standardmaessig. 5 erzeugt Untergruppe der Ordnung 2^{k-2} in (Z/2^k)*.

EN vs DE: Identisch.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine Fehler gefunden |

---

### =====================================================================
### PAPER 10 -- Wilson Quotient (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

1. **Basic values (Zeile 95-110 EN):** Nachgerechnet.
   - W_2 = (1!+1)/2 = 1. Korrekt.
   - W_3 = (2!+1)/3 = 1. Korrekt.
   - W_5 = (24+1)/5 = 5. Korrekt.
   - W_7 = (720+1)/7 = 103. Korrekt (7*103 = 721).
   - W_11 = (3628800+1)/11 = 329891. Korrekt (11*329891 = 3628801).

2. **Half-factorial formula (Zeilen 133-200 EN):** Korrekt.
   - (p-1)! = prod j(p-j) fuer j=1..m: Korrekt.
   - p-j = -j(1 - pj^{-1}): Korrekt.
   - prod(1-pj^{-1}) equiv 1 - p*sum(j^{-1}) mod p^2: Korrekt (Taylorentwicklung, hoehre Terme O(p^2)).
   - (-1)^m(m!)^2 equiv -1 mod p: Dies folgt aus Wilson. Korrekt.

3. **Beispiel p=5:** mu_5 = (1*4+1)/5 = 5/5 = 1. S_5 = 1+3 = 4 mod 5. W_5 equiv 1+4 = 5 equiv 0 mod 5. Korrekt.

4. **Beispiel p=7:** mu_7 = (-1*36+1)/7 = -35/7 = -5. S_7 = 1+4+5 = 10 equiv 3 mod 7. W_7 equiv -5+3 = -2 equiv 5 mod 7. Korrekt. (W_7 = 103 = 14*7 + 5, also 103 equiv 5 mod 7.)

5. **Bernoulli connection:** Beweisskizze. Verweist korrekt auf Ireland-Rosen und BEW1998.

| ID | Paper | Sprache | Schwere | Zeile | Beschreibung | Fix |
|----|-------|---------|---------|-------|-------------|-----|
| BUG-B2-P10-DE-001 | 10 | DE | GERING | 231 | Kommentar im LaTeX-Quelltext: "% BUG-B2-P10-DE-004: Schlussfolgerung korrigiert..." -- Interner Build-Kommentar, der im PDF nicht sichtbar ist, aber im Quelltext stehen bleibt. Unprofessionell fuer ein publiziertes Paper. | Kommentar entfernen |

---

### =====================================================================
### PAPER 11 -- Wilson Abelian Groups (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

1. **Pairing Lemma (Zeilen 86-102 EN):** Korrekt. Nicht-selbst-inverse heben sich paarweise auf.

2. **Subgroup Lemma (Zeilen 104-123 EN):** Korrekt. S = {g : g^2 = e} ist Untergruppe (Abschluss unter Multiplikation da abelsch, Inverse = Elemente selbst). S ist elementar-abelsch: jedes Element hat Ordnung 1 oder 2, also S = (Z/2Z)^k.

3. **General Wilson Theorem (Zeilen 129-179 EN):** Korrekt.
   - k=0: Produkt = e.
   - k=1: Produkt = tau.
   - k>=2: Alternativer Beweis via tau^{2^{k-1}} = e korrekt. Paarung s -> tau*s hat keine Fixpunkte (da tau != e). |S|/2 = 2^{k-1} Paare, jedes gibt tau. Produkt = tau^{2^{k-1}}, und 2^{k-1} gerade fuer k>=2, also tau^{2^{k-1}} = e. Korrekt.

4. **Gauss Generalisation (Zeilen 212-299 EN):** Korrekt. Analyse der Anzahl der Elemente der Ordnung <= 2 in (Z/nZ)* via CRT. Die Werte n = 1, 2, 4, p^k, 2p^k sind korrekt identifiziert.

5. **Numerical check n=8:** 1*3*5*7 = 105 = 13*8 + 1 equiv 1 mod 8. Korrekt.

EN vs DE: Identisch.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine Fehler gefunden |

---

### =====================================================================
### PAPER 12 -- Wilson Applications (EN + DE)
### =====================================================================

**Urteil: DRUCKREIF**

1. **Primality characterisation (Zeilen 79-110 EN):** Identisch mit Paper 8 Wilson-Beweis. Korrekt.

2. **-1 as QR (Zeilen 143-208 EN):** Korrekt.
   - (p-1)! = M^2 * (-1)^{(p-1)/2} wobei M = ((p-1)/2)!: Korrekt (pairing k <-> p-k, p-k equiv -k mod p).
   - Kombination mit Wilson: -1 equiv M^2 * (-1)^{(p-1)/2}.
   - p equiv 1 mod 4 => (-1)^{(p-1)/2} = 1, also -1 equiv M^2, -1 ist QR. Korrekt.
   - p equiv 3 mod 4 => (-1)^{(p-1)/2} = -1, also -1 equiv -M^2, also M^2 equiv 1. Falls -1 = x^2, dann (-1)^{(p-1)/2} = x^{p-1} = 1 (Fermat), aber (-1)^{(p-1)/2} = -1, Widerspruch. Korrekt.

3. **Wolstenholme harmonic (Zeilen 226-256 EN):** Korrekt. sum j^{-1} = sum j (da j -> j^{-1} ist Bijektion auf {1,...,p-1}). sum j = p(p-1)/2 equiv 0 mod p. Beachte: Dies ist das SCHWACHE Wolstenholme-Lemma (mod p), nicht das starke (mod p^2).

4. **Fermat via Wilson (Zeilen 279-319 EN):** Korrekt. Standard-Permutationsargument. Cancellation ist gueltig da (p-1)! nicht equiv 0 mod p.

5. **AKS-Verweis:** \cite{Agrawal2004} korrekt als "PRIMES is in P", Ann. Math. 160(2):781-793, 2004. Korrekt.

EN vs DE: Identisch.

| ID | Schwere | Beschreibung |
|----|---------|-------------|
| -- | -- | Keine Fehler gefunden |

---

## ZUSAMMENFASSUNG ALLER OFFENEN BUGS (Stand Build 14)

### HOCH
| ID | Paper | Beschreibung |
|----|-------|-------------|
| BUG-B1-P2-EN-001 | 2 EN | Fall p=5, q>=17 nur per Computersuche, kein elementarer Beweis |

### MITTEL
| ID | Paper | Beschreibung |
|----|-------|-------------|
| BUG-B1-P2-EN-002 | 2 EN | Schranke Fall p>=7 unzureichend begruendet |
| BUG-B1-P2-DE-001 | 2 DE | Algebraisch falscher Zwischenschritt in Lemma schranke_k, Zeile 217 |
| BUG-B14-P2-EN-004 | 2 EN | Kongruenz pq equiv alpha mod (p-1)(q-1) nicht explizit bewiesen |
| BUG-B14-P2-DE-002 | 2 DE | Falsche Begruendung fuer (p-1) | (q-1) in Zeile 250 |
| BUG-B14-P2-DE-003 | 2 DE | Lemma grenze_q nur fuer D=1 bewiesen, allgemeiner Fall fehlt |

### GERING
| ID | Paper | Beschreibung |
|----|-------|-------------|
| BUG-B2-P8-EN-001 | 8 EN+DE | "2k-2 >= k+1 > k" falsch fuer k=3 (Gleichheit statt >) |
| BUG-B2-P10-DE-001 | 10 DE | Interner Build-Kommentar im Quelltext |
| BUG-011 | 7 DE | Stilanmerkung: irrefuehrende Ungleichung in Remark |
| BUG-B14-P2-EN-003 | 2 EN | Minor: q >= p+2 Annahme koennte klarer sein |

### KOSMETISCH
| ID | Paper | Beschreibung |
|----|-------|-------------|
| -- | -- | Keine neuen kosmetischen Fehler |

---

## PAPER-UEBERSICHT

| Paper | Titel | EN | DE | Urteil |
|-------|-------|----|----|--------|
| 1 | Giuga 3-Prim | OK | OK | DRUCKREIF |
| 2 | Lehmer 3-Prim | 3 Bugs | 4 Bugs | UEBERARBEITUNG |
| 3 | Giuga-Carmichael | OK | OK | DRUCKREIF |
| 4 | Giuga quadratfrei | OK | OK | DRUCKREIF |
| 5 | Giuga kein Semiprim | OK | OK | DRUCKREIF |
| 6 | Lehmer quadratfrei | OK | OK | DRUCKREIF |
| 7 | Lehmer kein Semiprim | OK | 1 Stil | DRUCKREIF |
| 8 | Wilson Grundsatz | 1 gering | 1 gering | DRUCKREIF |
| 9 | Wilson Primzahlpot. | OK | OK | DRUCKREIF |
| 10 | Wilson Quotient | OK | 1 gering | DRUCKREIF |
| 11 | Wilson abelsch | OK | OK | DRUCKREIF |
| 12 | Wilson Anwendungen | OK | OK | DRUCKREIF |
