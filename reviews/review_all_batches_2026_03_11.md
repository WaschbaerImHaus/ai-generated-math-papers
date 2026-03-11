# Vollständiger Review-Bericht: Alle Batches (Papers 1–28)
**Datum:** 2026-03-11
**Build:** 7
**Gutachter:** Claude Code (Automatisierter mathematischer Review)
**Scope:** Formale Korrektheit, Schlüssigkeit, Vollständigkeit, EN↔DE-Vergleich

---

## ZUSAMMENFASSUNG

| Batch | Papers | EN-Urteil | DE-Urteil | Neue Bugs |
|---|---|---|---|---|
| Batch 1 | 1–7 (Giuga/Lehmer) | ✅ alle druckreif | ⚠️ 2 Stilanm. | BUG-010, 011 (bekannt) |
| Batch 2 | 8–12 (Wilson) | ✅ alle druckreif | ✅ alle druckreif | keine |
| Batch 3 | 13–16 (Sieblehre) | ✅ alle druckreif | ✅ alle druckreif | keine |
| Batch 4 | 17–20 (Kreismethode) | ⚠️ 3 EN-Fehler | ⚠️ 6 DE-Fehler | bekannte Bugs bestätigt |
| Batch 5 | 21–24 (Riemann/RH) | ✅ alle druckreif | ✅ alle druckreif | BUG-B5-P24 (Stil) |
| Batch 6 | 25–28 (BSD/Elliptisch) | ⚠️ 1 EN-Fehler | ❌ 6 DE-Fehler | BUG-B6-* (NEU) |

**Druckreife Papers gesamt:** 40/56
**Fehlerhafte Papers:** 16/56

---

## BATCH 1: Giuga & Lehmer (Papers 1–7)

### Paper 1 — Giuga: Keine 3-Primfaktor-Produkte als Giuga-Pseudoprimen
| | EN | DE |
|---|---|---|
| LaTeX-Struktur | ✅ | ✅ |
| Math. Korrektheit | ✅ | ✅ |
| Beweise vollständig | ✅ | ✅ |
| EN↔DE-Konsistenz | — | ✅ identisch |
| **Urteil** | **DRUCKREIF** | **DRUCKREIF** |

### Paper 2 — Lehmer: Totient-Problem für 3-Primfaktor-Produkte
| | EN | DE |
|---|---|---|
| LaTeX-Struktur | ✅ | ✅ |
| Math. Korrektheit | ✅ | ⚠️ |
| Beweise vollständig | ✅ | ✅ |
| EN↔DE-Konsistenz | — | ⚠️ Z.96–98 schwächer |
| **Urteil** | **DRUCKREIF** | **STILANMERKUNG** (BUG-010) |

**BUG-010:** Zeile 96–98 DE — Ungleichungsabschätzung schwächer als EN. Ergebnis korrekt, Formulierung unklar.

### Paper 3 — Giuga–Carmichael: Unverträglichkeit
- **EN + DE: DRUCKREIF.** CRT-Berechnung für {3,5,7} vollständig verifiziert.

### Paper 4 — Giuga: Quadratfrei
- **EN + DE: DRUCKREIF.** Kompakter, korrekter Beweis via Primteiler-Widerspruch.

### Paper 5 — Giuga: Kein Semiprim
- **EN + DE: DRUCKREIF.** Elegantes Ordnungsargument (p-1 < q → p-1 = 0 → p = 1, Widerspruch).

### Paper 6 — Lehmer: Quadratfrei
- **EN + DE: DRUCKREIF.** Beweis via φ(n)-Formel vollständig.

### Paper 7 — Lehmer: Kein Semiprim
| | EN | DE |
|---|---|---|
| Math. Korrektheit | ✅ | ⚠️ |
| Beweise vollständig | ✅ | ✅ |
| EN↔DE-Konsistenz | — | ⚠️ Z.96–98 Schranke falsch |
| **Urteil** | **DRUCKREIF** | **STILANMERKUNG** (BUG-011) |

**BUG-011:** Zeile 96–98 DE — Schranke schreibt `1·p-1` statt korrekt `1·(q-2)-1` oder `1·2-1`. Ergebnis stimmt, Herleitung nicht schlüssig.

---

## BATCH 2: Wilson-Theorem (Papers 8–12)

**Alle 10 Dateien: ✅ DRUCKREIF. Keine Fehler.**

- Paper 8 (Wilson-Theorem): Pairing-Argument elegant und korrekt.
- Paper 9 (Primzahlpotenzen): p=2-Sonderfall (4 selbst-inverse Elemente) korrekt behandelt.
- Paper 10 (Wilson-Quotient): Half-factorial formula vollständig. EN+DE konsistent.
- Paper 11 (Abelsche Gruppen): Gauss-Verallgemeinerung (n ∈ {1,2,4,p^k,2p^k}) vollständig bewiesen.
- Paper 12 (Anwendungen): 4 Anwendungen (Primtest, quad. Reste, Wolstenholme, Fermat) korrekt.

---

## BATCH 3: Sieblehre (Papers 13–16)

**Alle 8 Dateien: ✅ DRUCKREIF. Keine neuen Fehler.**

- Paper 13 (Brunns Satz): Bonferroni-Ungleichungen und Brun-Konstante korrekt.
- Paper 14 (Selbergs Sieb): λ²-Methode und Brun-Titchmarsh korrekt.
- Paper 15 (Große Sieb-Ungleichung): Beweisskizzen gekennzeichnet, mathematisch korrekt.
- Paper 16 (Chens Satz): Buchstab-Switching-Faktor 1/2 korrekt erklärt.

---

## BATCH 4: Kreismethode & Anwendungen (Papers 17–20)

### Paper 17 — Kreismethode (Hardy–Littlewood)
- **EN + DE: DRUCKREIF.** Farey-Zerlegung, Major/Minor Arcs vollständig. DE ist korrekte Übersetzung.

### Paper 18 — Vinogradovs Drei-Primzahlen-Satz
| Bug-ID | Sprache | Zeile | Beschreibung |
|---|---|---|---|
| BUG-B4-P18-EN-001 | EN | ~300 | Vorzeichen in Remark 4.3: c₂(n) via `e^{πin}` statt direkt. Für gerades n: c₂(n) = 1 nicht klar hergeleitet. |
| BUG-B4-P18-DE-001 | DE | ~203 | Identisch: `e^{πin} = 1` fehlerhaft für gerades n. |
| BUG-B4-P18-DE-002 | DE | ~226 | Beweis von Satz 3.3 (singuläres Integral) stoppt nach Simplex-Diskussion. Fallunterscheidung unvollständig skizziert. |

**Urteile:** EN: FEHLER | DE: FEHLER

### Paper 19 — Waringsches Problem
| Bug-ID | Sprache | Zeile | Beschreibung |
|---|---|---|---|
| BUG-B4-P19-EN-001 | EN | ~208 | Lemma 5.5 aus Paper 17 falsch zitiert: J(n) als Konstante dargestellt, korrekt wäre J(n) ~ C·n^{s/k-1}. |
| BUG-B4-P19-DE-001 | DE | ~161 | Identischer Fehler wie EN-001. |
| BUG-B4-P19-DE-002 | DE | ~244 | Definition "p-adische Kongruenzrechnung" fehlt (vor Theorem 3.1). EN hat sie (Z.305-309). |
| BUG-B4-P19-DE-003 | DE | ~209 | Beweis von Satz 2.3 ist nur Skizze; zentraler Schritt nicht begründet. |

**Urteile:** EN: FEHLER | DE: FEHLER

### Paper 20 — Goldbach'sche singuläre Reihe
| Bug-ID | Sprache | Zeile | Beschreibung |
|---|---|---|---|
| BUG-B4-P20-EN-001 | EN | ~501 | Satzfragment in Remark (Parity Barrier): "...and the analogue for even n" bricht ab. |
| BUG-B4-P20-EN-002 | EN | ~246 | Tabelle zeigt S(n) statt vorhergesagtes r₂(n) = S(n)·n/(log n)². Spalte "Prognose" fehlerhaft. |
| BUG-B4-P20-DE-001 | DE | ~185 | Identische Tabellenfehler wie EN-002. |

**Urteile:** EN: FEHLER | DE: FEHLER

---

## BATCH 5: Riemann-Zeta & RH-Ansätze (Papers 21–24)

**Dieser Batch wurde in diesem Build erstmalig begutachtet.**

### Paper 21 — Riemann-Zeta-Funktion
- **EN + DE: ✅ DRUCKREIF.**
- Euler-Produkt, Funktionalgleichung, analytische Fortsetzung: alles korrekt.
- EN↔DE vollständig konsistent.

### Paper 22 — Nichttriviale Nullstellen
- **EN + DE: ✅ DRUCKREIF.**
- Riemann-von-Mangoldt-Formel N(T) korrekt. γ₁ = 14.134725... korrekt.
- GUE-Paarkorrelation r₂(s) = 1 - (sin(πs)/πs)² korrekt.

### Paper 23 — Explizite Formel und Primzahlsatz
- **EN + DE: ✅ DRUCKREIF.**
- Chebyshev-Funktionen, Perron-Formel, Schoenfeld-Schranke korrekt.
- Explizite Formel für ψ(x) vollständig.

### Paper 24 — RH-Ansätze
| Bug-ID | Sprache | Zeile | Beschreibung |
|---|---|---|---|
| BUG-B5-P24-EN-001 | EN | 38–39, 139 | Inkonsistente Formatierung "Hilbert–Pólya" (Gedankenstrich statt --). Nur kosmetisch. |
| BUG-B5-P24-DE-001 | DE | 38–39, 139 | Identisch. |

**Urteile:** EN: DRUCKREIF (Stilanm.) | DE: DRUCKREIF (Stilanm.)

Mathematisch: Berry-Keating-Hamiltonian, GUE-Statistik, Selberg-Spurformel, Levinson-Conrey (0.41) alle korrekt.

---

## BATCH 6: Elliptische Kurven & BSD (Papers 25–28)

**Dieser Batch wurde in diesem Build erstmalig begutachtet.**

### Paper 25 — Elliptische Kurven über Q
| | EN | DE |
|---|---|---|
| Gruppengesetz | ✅ Chord-and-tangent korrekt | ✅ |
| Mazur-Torsionssatz (15 Gruppen) | ✅ | ✅ |
| Mordell-Weil | ✅ | ✅ |
| Néron-Tate Höhe | ✅ | ⚠️ |
| **Urteil** | **DRUCKREIF** | **FEHLER** |

**BUG-B6-P25-DE-001** (Minor): DE fehlt das Beispiel für singuläre Kurve y²=x³-x² (EN Zeile 231–232). DE springt direkt zu y²=x³+1.

### Paper 26 — L-Funktion elliptischer Kurven
| | EN | DE |
|---|---|---|
| Frobenius-Spur a_p | ✅ | ✅ |
| Euler-Produkt (3 Fälle) | ✅ | ✅ |
| Verifikationstabelle a_5 für E:y²=x³-x | ✅ vollständig | ❌ fehlt |
| Modularitätssatz (Wiles 1995) | ✅ mit Fermat-Korollar | ⚠️ ohne Fermat-Erwähnung |
| **Urteil** | **DRUCKREIF** | **FEHLER** |

**BUG-B6-P26-DE-001** (Major): Verifikationstabelle x→x³-x mod 5 fehlt vollständig in DE (EN Zeile 120–128).
**BUG-B6-P26-DE-002** (Minor): DE erwähnt nicht, dass Fermats Letzter Satz als Korollar des Modularitätssatzes folgt (via Ribet).

### Paper 27 — BSD-Vermutung
| | EN | DE |
|---|---|---|
| Schwache BSD (Kolyvagin 1990) | ✅ | ✅ |
| Starke BSD-Formel | ✅ alle Ingredienzen | ⚠️ |
| Tate-Shafarevich-Gruppe | ✅ | ✅ |
| Coates-Wiles (1977) | ✅ | ✅ |
| **Urteil** | **DRUCKREIF** | **FEHLER** |

**BUG-B6-P27-DE-001** (Minor): DE fehlt Aussage "Strong BSD wurde numerisch für viele Kurven verifiziert, aber nur in Spezialfällen bewiesen" (EN Zeile 162–173).

### Paper 28 — Kongruente Zahlen und BSD
| | EN | DE |
|---|---|---|
| Definition kongruente Zahl | ✅ | ✅ |
| Äquivalenz E_n: y²=x³-n²x | ✅ Parametrisierung korrekt | ❌ fehlerhafte Param. |
| Tunnell-Theorem (1983) | ✅ A(n),B(n),C(n),D(n) korrekt | ❌ C(n),D(n) falsch |
| Torsion Z/2Z×Z/2Z | ⚠️ Nagell-Lutz unklar | ✅ |
| Verifikation n=5,6,7 | ✅ | ✅ |
| **Urteil** | **FEHLER** | **KRITISCHER FEHLER** |

**BUG-B6-P28-EN-001** (Major): Nagell-Lutz-Bedingung falsch formuliert. Für E_n: a=-n², b=0 ist 4a³+27b²=-4n⁶ negativ; Formulierung "y₀²|-4n⁶" ist technisch inkorrekt. Korrekt: y₀² | 4n⁶.

**BUG-B6-P28-DE-001** (Kritisch): Beweis Theorem 2.1 verwendet fehlerhafte Parametrisierung: x₀=((a+b)/2)² statt korrekt x₀=c²/4. Text enthält sogar ein eingebettetes "nein --" als Hinweis auf Selbstkorrektur-Versuch, ohne dies abzuschließen.

**BUG-B6-P28-DE-002** (Minor): Formulierung "achtes Pythagoraisches Tripel" ist unverständlich. EN sagt klar "a ≠ b (otherwise the triangle degenerates)".

**BUG-B6-P28-DE-003** (KRITISCH — MATHEMATISCHER FEHLER):
Tunnell-Formeln für C(n) und D(n) haben falsche Koeffizienten:
- EN (korrekt): `C(n) = #{4x²+y²+8z²=n/2}`, `D(n) = #{4x²+y²+32z²=n/2}`
- DE (falsch):   `C(n) = #{8x²+2y²+16z²=n/2}`, `D(n) = #{8x²+2y²+64z²=n/2}`

Die DE-Variante ist eine falsch skalierte Version (×2 in allen Koeffizienten), die mathematisch nicht dem Tunnell-Kriterium entspricht.

---

## GESAMTÜBERSICHT ALLER OFFENEN BUGS

### Kritische / Hohe Schwere

| Bug-ID | Paper | Sprache | Art | Beschreibung |
|---|---|---|---|---|
| BUG-B4-P20-EN-002 | 20 | EN | Rechenfehler | Tabelle S(n) systematisch falsch |
| BUG-B4-P20-DE-001 | 20 | DE | Rechenfehler | Identisch EN-002 |
| BUG-B6-P28-DE-001 | 28 | DE | Math. Fehler | Fehlerhafte Parametrisierung im Beweis |
| BUG-B6-P28-DE-003 | 28 | DE | Math. Fehler | Tunnell-Formel Koeffizienten falsch (×2-Skalierung) |

### Mittlere Schwere

| Bug-ID | Paper | Sprache | Art | Beschreibung |
|---|---|---|---|---|
| BUG-B4-P18-EN-001 | 18 | EN | Math. Inkonsistenz | Remark 4.3: c₂(n) via e^{πin} unklar |
| BUG-B4-P18-DE-001 | 18 | DE | Math. Inkonsistenz | Identisch EN-001 |
| BUG-B4-P18-DE-002 | 18 | DE | Unvollst. Beweis | Singuläres Integral Fallunterscheidung fehlt |
| BUG-B4-P19-EN-001 | 19 | EN | Zitierfehler | J(n) als Konstante statt ~ C·n^{s/k-1} |
| BUG-B4-P19-DE-001 | 19 | DE | Zitierfehler | Identisch EN-001 |
| BUG-B4-P19-DE-002 | 19 | DE | Fehlende Def. | p-adische Kongruenzrechnung nicht definiert |
| BUG-B4-P19-DE-003 | 19 | DE | Unvollst. Beweis | Satz 2.3 nur Skizze |
| BUG-B6-P26-DE-001 | 26 | DE | Fehlender Inhalt | Verifikationstabelle a_5 fehlt |
| BUG-B6-P28-EN-001 | 28 | EN | Formelfehler | Nagell-Lutz-Bedingung (y₀²|-4n⁶ statt 4n⁶) |

### Geringe Schwere / Stilistisch

| Bug-ID | Paper | Sprache | Art | Beschreibung |
|---|---|---|---|---|
| BUG-010 | 2 | DE | Stil | Strukturelle Redundanz (2 Beweise ohne Trennung) |
| BUG-011 | 7 | DE | Stil | Schwacher Zwischenschritt (1·p-1 statt 1·(q-2)-1) |
| BUG-B4-P20-EN-001 | 20 | EN | Satzfragment | "...and the analogue for even n" bricht ab |
| BUG-B5-P24-EN-001 | 24 | EN | Stil | Hilbert-Pólya Bindestrich-Formatierung |
| BUG-B5-P24-DE-001 | 24 | DE | Stil | Identisch EN-001 |
| BUG-B6-P25-DE-001 | 25 | DE | Fehlend | Beispiel singuläre Kurve fehlt |
| BUG-B6-P26-DE-002 | 26 | DE | Fehlend | Fermat-Korollar des Modularitätssatzes nicht erwähnt |
| BUG-B6-P27-DE-001 | 27 | DE | Fehlend | Fehlende Aussage zur numerischen Verifikation |
| BUG-B6-P28-DE-002 | 28 | DE | Formulierung | "achtes Pythagoraisches Tripel" unklar |

---

## HANDLUNGSEMPFEHLUNGEN (nach Priorität)

### Sofort-Korrekturen (vor nächstem Build)
1. **BUG-B6-P28-DE-003:** Tunnell-Formel C(n), D(n) auf korrekte Koeffizienten setzen.
2. **BUG-B6-P28-DE-001:** Parametrisierung x₀=c²/4 aus EN übernehmen, "nein --"-Fragment entfernen.
3. **BUG-B4-P20-EN-002 + DE-001:** Tabelle S(n) korrekt neu berechnen oder Spalte als "Werte S(n)" statt "Prognose r₂(n)" kennzeichnen.

### Vor Druckfreigabe
4. **BUG-B4-P18-EN-001 + DE-001:** Remark 4.3 — c₂(n) direkt definieren ohne e^{πin}.
5. **BUG-B4-P18-DE-002:** Fallunterscheidung für singuläres Integral aus EN ergänzen.
6. **BUG-B4-P19-EN-001 + DE-001:** J(n) ~ C·n^{s/k-1} korrekt formulieren.
7. **BUG-B4-P19-DE-002:** Definition p-adische Kongruenzrechnung vor Theorem 3.1 einfügen.
8. **BUG-B6-P26-DE-001:** Verifikationstabelle aus EN übernehmen.
9. **BUG-B6-P28-EN-001:** Nagell-Lutz-Bedingung auf y₀² | 4n⁶ korrigieren.

### Optional / Stilistisch
10. BUG-010, BUG-011: Kleine Stilverbesserungen in Batch 1 DE.
11. BUG-B4-P20-EN-001: Satzfragment vervollständigen.
12. BUG-B5-P24-EN/DE-001: Hilbert-Pólya Schreibweise harmonisieren.
13. BUG-B6-P25-DE-001: Beispiel singuläre Kurve ergänzen.
14. BUG-B6-P26-DE-002: Fermat-Korollar erwähnen.
15. BUG-B6-P27-DE-001: Hinweis auf numerische Verifikation ergänzen.
16. BUG-B6-P28-DE-002: "achtes Pythagoraisches Tripel" ersetzen.
