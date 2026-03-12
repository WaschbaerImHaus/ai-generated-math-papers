# bruns_constant.py — Bruns Konstante B₂ und Zwillingsprimzahlen

**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-12
**Modul:** `src/bruns_constant.py`

---

## Inhaltsverzeichnis

1. [Mathematischer Hintergrund](#mathematischer-hintergrund)
2. [Klassen- und Methodenübersicht](#klassen--und-methodenübersicht)
3. [Wichtigste Algorithmen](#wichtigste-algorithmen)
4. [Beispielanwendungen](#beispielanwendungen)

---

## Mathematischer Hintergrund

### Zwillingsprimzahlen

Ein **Zwillingsprimzahlpaar** ist ein Paar $(p, p+2)$, bei dem beide Zahlen prim sind:

$$(3, 5),\ (5, 7),\ (11, 13),\ (17, 19),\ (29, 31),\ (41, 43), \ldots$$

**Zwillingsprim-Vermutung** (offen): Es gibt unendlich viele Zwillingsprimzahlpaare.

Die Frage nach der Unendlichkeit dieser Paare ist eines der ältesten offenen Probleme der Zahlentheorie.

### Bruns Konstante

**Viggo Brun** bewies 1919 einen bemerkenswerten Satz:

> **Bruns Theorem (1919):** Die Reihe der reziproken Zwillingsprimzahlen konvergiert:
> $$B_2 = \sum_{\substack{(p,\, p+2) \\ \text{Zwillingsprim}}} \left(\frac{1}{p} + \frac{1}{p+2}\right) < \infty$$

Dies ist bemerkenswert, denn die harmonische Reihe aller Primzahlreziproken divergiert ($\sum 1/p = \infty$), aber die Teilreihe über Zwillingsprimzahlen konvergiert — unabhängig davon, ob es unendlich oder endlich viele Zwillingsprimzahlpaare gibt.

**Numerischer Wert:**
$$B_2 \approx 1.9021605831040309591714\ldots$$

Berechnet bis $p < 10^{14}$ (Nicely 1995, Sebah & Gourdon 2002).

### Partielle Bruns Summe

$$B_2(x) = \sum_{\substack{(p, p+2) \text{ Zwillingsprim} \\ p \leq x}} \left(\frac{1}{p} + \frac{1}{p+2}\right)$$

$B_2(x)$ konvergiert sehr langsam gegen $B_2$. Der Fehler beträgt grob:

$$B_2 - B_2(x) = O\!\left(\frac{1}{(\ln x)^2}\right)$$

### Hardy-Littlewood Conjecture B

**Hardy & Littlewood (1923)** vermuteten für die Anzahl $\pi_2(x)$ der Zwillingsprimzahlpaare bis $x$:

$$\pi_2(x) \sim 2 C_2 \cdot \frac{x}{(\ln x)^2} \quad (x \to \infty)$$

mit der **Zwillingsprimzahlkonstante**:
$$C_2 = \prod_{\substack{p > 2 \\ p \text{ prim}}} \frac{p(p-2)}{(p-1)^2} \approx 0.6601618158\ldots$$

**Status:** Unbewiesen (Vermutung). Die ersten Faktoren des Produkts sind:

| $p$ | Faktor $\frac{p(p-2)}{(p-1)^2}$ | Teilprodukt |
|-----|----------------------------------|-------------|
| 3 | $3/4 = 0.75$ | 0.75 |
| 5 | $15/16 = 0.9375$ | 0.703125 |
| 7 | $35/36 \approx 0.9722$ | 0.683594 |
| 11 | $99/100 = 0.99$ | 0.676758 |

**Bekannte exakte Werte:**

| $x$ | $\pi_2(x)$ |
|-----|------------|
| $10^4$ | 205 |
| $10^5$ | 1224 |
| $10^6$ | 8169 |
| $10^8$ | 440312 |
| $10^{10}$ | 27412679 |

---

## Klassen- und Methodenübersicht

### Klasse `BrunsKonstante`

Klassenattribut:
- `HARDY_LITTLEWOOD_C2 = 0.6601618158...` — die Zwillingsprimzahlkonstante $C_2$

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `berechne_zwillingsprimes` | `(grenze: int) → List[Tuple[int,int]]` | Alle Zwillingsprimzahlpaare $(p, p+2)$ mit $p \leq$ grenze |
| `bruns_summe` | `(grenze: int) → float` | Partielle Bruns Summe $B_2(x)$ |
| `konvergenzanalyse` | `(schranken_list?) → Dict` | Konvergenzprozess für verschiedene Schranken |
| `hochpraezise_berechnung` | `(grenze?, dezimalstellen?) → Dict` | Hochpräzise Berechnung via mpmath |
| `meissel_lehmer_abschaetzung` | `(x: float) → Dict` | Hardy-Littlewood-Abschätzung $\pi_2(x)$ |
| `hardylittlewood_vorhersage` | `(x: float) → Dict` | $C_2$-Berechnung und $\pi_2(x)$-Vorhersage |
| `vergleich_bruns_konstante_approximation` | `() → Dict` | Vergleich verschiedener Approximationsmethoden |

---

## Wichtigste Algorithmen

### 1. Sieb des Eratosthenes für Zwillingsprimzahlen

```
Eingabe: grenze (obere Schranke für kleineres Element p)
1. Erstelle bytearray ist_prim[0..grenze+2], initialisiert mit 1
2. Setze ist_prim[0] = ist_prim[1] = 0
3. Für i = 2, ..., floor(sqrt(grenze+2)):
       falls ist_prim[i]: markiere alle Vielfachen i²,i²+i,... als 0
4. Sammle alle p mit ist_prim[p]=1 und ist_prim[p+2]=1
```

**Komplexität:** $O(n \log \log n)$ für das Sieb.

### 2. Partielle Bruns Summe (einfache Summation)

```python
summe = 0.0
for p, p2 in berechne_zwillingsprimes(grenze):
    summe += 1.0/p + 1.0/p2
return summe
```

### 3. Hochpräzise Summation via mpmath

Für höhere Präzision wird mpmath mit `mp.dps = dezimalstellen + 10` verwendet, und jeder Term $1/p + 1/(p+2)$ wird als `mpf`-Zahl berechnet.

### 4. Richardson-Extrapolation

Aus zwei Summenwerten $B_2(x_1)$ und $B_2(x_2)$ wird ein besserer Schätzwert für $B_2$ ermittelt. Angesetzt wird das asymptotische Modell $B_2(x) \approx B_2 - C/\ln(x)$:

$$B_2 \approx \frac{B_2(x_2) \cdot \ln x_2 - B_2(x_1) \cdot \ln x_1}{\ln x_2 - \ln x_1}$$

Mit $x_1 = 50000$ und $x_2 = 100000$ liefert dies eine deutlich bessere Approximation als die direkte Summation.

### 5. Asymptotische Korrektur

Der Fehler bei Abschneiden bei $x$ wird abgeschätzt durch:

$$B_2 - B_2(x) \approx 4 C_2 \cdot \int_x^{\infty} \frac{dt}{t \cdot (\ln t)^2} \approx \frac{4 C_2}{\ln x}$$

Diese asymptotische Korrektur verbessert die Approximation erheblich.

### 6. Berechnung von $C_2$ als Teilprodukt

$$C_2 = \prod_{p=3}^{1000} \frac{p(p-2)}{(p-1)^2}$$

Das Produkt konvergiert schnell (jeder Faktor ist nahe 1), und der Fehler durch Abschneiden bei $p < 1000$ ist vernachlässigbar.

---

## Beispielanwendungen

### Erste Zwillingsprimpaare ausgeben

```python
from bruns_constant import BrunsKonstante

bk = BrunsKonstante()
paare = bk.berechne_zwillingsprimes(100)
# [(3,5), (5,7), (11,13), (17,19), (29,31), (41,43), (59,61), (71,73)]
print(f"Anzahl bis 100: {len(paare)}")
```

### Bruns Summe berechnen

```python
b2_1000 = bk.bruns_summe(1000)
b2_10000 = bk.bruns_summe(10000)
print(f"B₂(10³)  = {b2_1000:.10f}")
print(f"B₂(10⁴)  = {b2_10000:.10f}")
print(f"Bekannt: 1.9021605831")
```

### Konvergenzanalyse

```python
ka = bk.konvergenzanalyse([100, 1_000, 10_000, 100_000])
for r in ka['ergebnisse']:
    print(f"x={r['x']:8d}: B₂={r['b2_x']:.10f}, Paare={r['anzahl_paare']}")
```

### Hochpräzise Berechnung

```python
hp = bk.hochpraezise_berechnung(grenze=1_000_000, dezimalstellen=30)
print(f"B₂(10⁶) = {hp['b2_hochpräzise']}")
print(f"Laufzeit: {hp['laufzeit_sek']} Sek.")
```

### Hardy-Littlewood-Vorhersage

```python
hl = bk.hardylittlewood_vorhersage(1e8)
print(f"C₂ = {hl['c2_berechnet']:.10f}")   # ≈ 0.6601618158
print(f"π₂(10⁸) ≈ {hl['vorhersage_pi2_x']:.0f}")  # ≈ 440000

# Bekannte π₂-Werte zum Vergleich:
for exp in [4, 6, 8]:
    ml = bk.meissel_lehmer_abschaetzung(10**exp)
    print(f"π₂(10^{exp}) ≈ {ml['hl_abschaetzung']:.0f}")
```

### Approximationsvergleich

```python
approx = bk.vergleich_bruns_konstante_approximation()
print(f"Direkte Summation (10⁵): {approx['direkte_summation'][100_000]:.10f}")
print(f"Richardson-Extrapolation: {approx['richardson_extrapolation']['wert']:.10f}")
print(f"Asympt. Korrektur:        {approx['asymptotische_korrektur']['b2_approx']:.10f}")
print(f"Bekannt:                  {approx['bekannter_wert']:.10f}")
```

---

## Historische Anmerkungen

- **Viggo Brun (1919):** Beweis der Konvergenz von $\sum 1/p$ über Zwillingsprimzahlen via Bruns Sieb. Das Sieb liefert $\pi_2(x) = O(x (\log \log x)^2 / (\ln x)^2)$.
- **G.H. Hardy & J.E. Littlewood (1923):** Conjecture B, die asymptotische Formel $\pi_2(x) \sim 2C_2 x/(\ln x)^2$.
- **Thomas R. Nicely (1995):** Berechnete $B_2$ bis $p < 10^{14}$ — und entdeckte dabei den Pentium-FDIV-Bug!
- **Pascal Sebah & Xavier Gourdon (2002):** Verfeinerte Berechnung auf 12 gesicherte Dezimalstellen.

---

## Konvergenzbeweis (Idee)

Bruns Sieb liefert:
$$\pi_2(x) = O\!\left(\frac{x (\log \log x)^2}{(\ln x)^2}\right)$$

Daraus folgt für die Teilsummen (Abel-Summation):
$$B_2(x) = \int_3^x \frac{d\pi_2(t)}{t} + \int_3^x \frac{d\pi_2(t)}{t+2}$$

Die Konvergenz folgt, weil $\pi_2(t)/t \to 0$ schnell genug. Im Gegensatz: Bei gewöhnlichen Primzahlen gilt $\pi(t)/t \sim 1/\ln t$, und $\sum 1/p$ divergiert (Mertens).

---

*Dokumentation generiert für das specialist-maths Projekt. Autor: Michael Fuhrmann.*
