# Elliptische Kurven – Vollständige Dokumentation

**Datei:** `elliptic_curves.py`
**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10

---

## Inhaltsverzeichnis

1. [Grundlagen der elliptischen Kurven](#1-grundlagen)
2. [Gruppenstruktur und Additionsgesetz](#2-gruppenstruktur)
3. [Diskriminante und j-Invariante](#3-diskriminante-und-j-invariante)
4. [Implementierte Klassen](#4-implementierte-klassen)
5. [Zahlentheoretische Anwendungen](#5-zahlentheoretische-anwendungen)
6. [BSD-Vermutung](#6-bsd-vermutung)
7. [Shimura-Taniyama-Wiles und Fermats Letzter Satz](#7-shimura-taniyama-wiles)
8. [Kryptographie (ECC, ECDH)](#8-elliptic-curve-cryptography)
9. [Lenstra's ECM-Faktorisierung](#9-ecm-faktorisierung)
10. [Kongruente Zahlen](#10-kongruente-zahlen)
11. [Bekannte Kurven](#11-bekannte-kurven)

---

## 1. Grundlagen

Eine **elliptische Kurve** über einem Körper $K$ ist eine glatte, projektive algebraische Kurve vom Geschlecht 1 mit einem ausgezeichneten Punkt.

In der **Weierstraß-Normalform** schreibt man:

$$E: y^2 = x^3 + ax + b, \quad a, b \in K$$

Die Kurve ist **nicht-singulär** (also eine echte elliptische Kurve), wenn die **Diskriminante** von Null verschieden ist:

$$\Delta = -16(4a^3 + 27b^2) \neq 0$$

Singularität tritt auf wenn:
- $\Delta = 0$: Die Kurve hat einen **Knoten** ($4a^3 + 27b^2 = 0$) oder eine **Spitze**

### Beispiele

| Kurve | a | b | $\Delta$ | Bemerkung |
|-------|---|---|----------|-----------|
| $y^2 = x^3 - x$ | $-1$ | $0$ | $64$ | BSD-Kurve, 2-Torsion |
| $y^2 = x^3 + 1$ | $0$ | $1$ | $-1728$ | CM durch $\mathbb{Z}[\omega]$ |
| $y^2 = x^3 + 7$ | $0$ | $7$ | $-21952$ | secp256k1 (Bitcoin) |
| $y^2 = x^3 - n^2x$ | $-n^2$ | $0$ | — | Kongruente Zahl Kurve |

---

## 2. Gruppenstruktur

### Abelsche Gruppe

Die Menge der Punkte $E(K)$ zusammen mit dem **Punkt im Unendlichen** $\mathcal{O}$ bildet eine **abelsche Gruppe** $\langle E(K), +, \mathcal{O} \rangle$:

- **Neutralelement:** $\mathcal{O}$ (projektiver Punkt im Unendlichen)
- **Inverses:** $-P = (x, -y)$ für $P = (x, y)$ (Spiegelung an der $x$-Achse)
- **Addition:** Chord-and-Tangent-Methode

### Chord-and-Tangent-Methode

Die geometrische Interpretation:

1. Verbinde $P$ und $Q$ durch eine Gerade (Chord)
2. Diese trifft die Kurve in einem dritten Punkt $R'$
3. Spiegle $R'$ an der $x$-Achse: $R = -R' = P + Q$

Bei $P = Q$ verwende die Tangente an $P$ statt der Sekante.

### Additionsformeln

**Fall 1:** $P = \mathcal{O}$:
$$P + Q = Q$$

**Fall 2:** $Q = \mathcal{O}$:
$$P + Q = P$$

**Fall 3:** $P = -Q$ (inverse Punkte, $x_1 = x_2$, $y_1 = -y_2$):
$$P + Q = \mathcal{O}$$

**Fall 4:** $P = Q$ (Tangente, Punktverdopplung):
$$\lambda = \frac{3x_1^2 + a}{2y_1}$$

**Fall 5:** $P \neq Q$ (Sekante):
$$\lambda = \frac{y_2 - y_1}{x_2 - x_1}$$

**In beiden Fällen (4 und 5):**
$$x_3 = \lambda^2 - x_1 - x_2$$
$$y_3 = \lambda(x_1 - x_3) - y_1$$

### Skalarmultiplikation

Das $n$-fache eines Punktes $P$ wird durch den **Double-and-Add-Algorithmus** berechnet:

$$nP = \underbrace{P + P + \cdots + P}_{n \text{ mal}}$$

Durch Binärdarstellung von $n$ wird dies in $\mathcal{O}(\log n)$ Schritten berechnet:

```
R ← O
for bit in bits(n) from MSB to LSB:
    R ← 2R        (doppeln)
    if bit == 1:
        R ← R + P  (addieren)
return R
```

---

## 3. Diskriminante und j-Invariante

### Diskriminante

$$\Delta = -16(4a^3 + 27b^2)$$

Das Vorzeichen und der Wert bestimmen das qualitative Verhalten der Kurve:
- $\Delta > 0$: Die Kurve hat zwei Komponenten (ein ovales Stück und eine unbeschränkte Komponente)
- $\Delta < 0$: Die Kurve ist zusammenhängend (eine unbeschränkte Komponente)

### j-Invariante

$$j(E) = -1728 \cdot \frac{(4a)^3}{\Delta} = \frac{-1728 \cdot 64 a^3}{-16(4a^3 + 27b^2)} = \frac{1728 \cdot 4a^3}{4a^3 + 27b^2}$$

Die $j$-Invariante klassifiziert elliptische Kurven über dem algebraischen Abschluss:

$$E_1 \cong E_2 \text{ über } \bar{K} \iff j(E_1) = j(E_2)$$

Spezialwerte:
- $j = 0$: Kurve hat Automorphismengruppe der Ordnung 6 (CM durch $\mathbb{Z}[\omega]$, $\omega = e^{2\pi i/3}$)
- $j = 1728$: Kurve hat Automorphismengruppe der Ordnung 4 (CM durch $\mathbb{Z}[i]$)

---

## 4. Implementierte Klassen

### `ECPoint`

Repräsentiert einen Punkt $P = (x, y)$ auf einer elliptischen Kurve oder den Punkt im Unendlichen $\mathcal{O}$.

```python
# Verwendung
curve = EllipticCurve(-1, 0)        # y² = x³ - x
P = ECPoint(0, 0, curve)             # Punkt (0,0)
O = ECPoint.infinity(curve)          # Punkt im Unendlichen

Q = P + P                            # Verdopplung → O (da y=0)
R = 3 * P                            # Skalarmultiplikation
```

**Methoden:**

| Methode | Beschreibung |
|---------|-------------|
| `ECPoint.infinity(curve)` | Erstellt $\mathcal{O}$ |
| `is_infinity` | Prüft ob $P = \mathcal{O}$ |
| `__add__(Q)` | Gruppenaddition $P + Q$ |
| `__neg__()` | Inverses $-P = (x, -y)$ |
| `__mul__(n)` | Skalarmultiplikation $nP$ |
| `on_curve()` | Prüft ob $P \in E(K)$ |

### `EllipticCurve`

Elliptische Kurve über $\mathbb{R}$.

```python
curve = EllipticCurve(-1, 0)
print(curve.j_invariant())           # 1728.0
print(curve.discriminant())          # 64.0
print(curve.is_on_curve(1, 0))      # True
points = curve.points_over_fp(7)    # E(F_7)
```

### `EllipticCurveModP`

Elliptische Kurve über $\mathbb{F}_p = \mathbb{Z}/p\mathbb{Z}$.

```python
curve = EllipticCurveModP(-1, 0, 7)
points = curve.all_points()          # Alle Punkte in E(F_7)
n = curve.group_order()              # #E(F_7)
```

### `ECCKeyExchange`

ECDH-Schlüsselaustausch.

```python
curve = EllipticCurveModP(-1, 0, 101)
G = curve.all_points()[1]            # Basispunkt
ecdh = ECCKeyExchange(curve, G)
result = ecdh.demo_exchange()        # Vollständige Demo
```

---

## 5. Zahlentheoretische Anwendungen

### Satz von Hasse

Für jede Primzahl $p$ und elliptische Kurve $E/\mathbb{F}_p$ gilt:

$$|{#}E(\mathbb{F}_p) - (p+1)| \leq 2\sqrt{p}$$

Dies bedeutet: $\#E(\mathbb{F}_p) \approx p+1$ mit einem Fehler von höchstens $2\sqrt{p}$.

### Frobenius-Spur

$$a_p = p + 1 - \#E(\mathbb{F}_p)$$

Die Frobenius-Spur ist der zentrale Parameter der $L$-Funktion von $E$.

### Nagell-Lutz-Satz

Für eine elliptische Kurve $E: y^2 = x^3 + ax + b$ mit $a, b \in \mathbb{Z}$ gilt:

**Satz** (Nagell 1935, Lutz 1937): Ist $P = (x, y)$ ein **Torsionspunkt** in $E(\mathbb{Q})$, dann:
1. $x, y \in \mathbb{Z}$ (ganzzahlige Koordinaten)
2. $y = 0$ oder $y^2 \mid \Delta$ (y² teilt die Diskriminante)

---

## 6. BSD-Vermutung

Die **Birch und Swinnerton-Dyer-Vermutung** (BSD) ist eines der sieben Millennium-Probleme und verbindet die analytische mit der algebraischen Theorie elliptischer Kurven.

### Mordell-Weil-Gruppe

Nach dem **Satz von Mordell-Weil**:

$$E(\mathbb{Q}) \cong \mathbb{Z}^r \oplus E(\mathbb{Q})_{\text{tors}}$$

wobei:
- $r \geq 0$: **Mordell-Weil-Rang** (Anzahl unabhängiger rationaler Punkte unendlicher Ordnung)
- $E(\mathbb{Q})_{\text{tors}}$: Endliche **Torsionsgruppe** (nach Mazur: max. 16 Elemente)

### L-Funktion

Die L-Funktion von $E$ wird definiert durch das **Euler-Produkt**:

$$L(E, s) = \prod_{p \nmid \Delta} \frac{1}{1 - a_p p^{-s} + p^{1-2s}} \cdot \prod_{p \mid \Delta} (\text{schlechte Reduktion})$$

### BSD-Vermutung (grob)

$$\text{ord}_{s=1} L(E, s) = r = \text{Rang}(E(\mathbb{Q}))$$

Das heißt: **Die Ordnung der Nullstelle von $L(E,s)$ bei $s=1$ ist gleich dem Rang.**

Weiterhin macht BSD eine präzise Aussage über den Leitkoeffizienten:

$$\lim_{s \to 1} \frac{L(E,s)}{(s-1)^r} = \frac{\Omega \cdot \text{Reg} \cdot \prod_p c_p \cdot |\text{III}|}{|E(\mathbb{Q})_{\text{tors}}|^2}$$

---

## 7. Shimura-Taniyama-Wiles

### Der Satz

**Shimura-Taniyama-Weil-Vermutung** (bewiesen von Wiles 1995, Taylor-Wiles):

> Jede semistabile elliptische Kurve über $\mathbb{Q}$ ist **modular**, d.h. sie ist ein Quotient der Jacobi-Varietät einer Modulkurve.

Konkret: Für jede elliptische Kurve $E/\mathbb{Q}$ existiert eine Modulform $f$ vom Gewicht 2, Stufe $N$ (dem Führer von $E$), sodass:

$$L(E, s) = L(f, s)$$

### Verbindung zu Fermats Letztem Satz

Andrew Wiles' Beweis des **Großen Fermatschen Satzes** (1995) nutzt diese Verbindung:

1. **Frey** (1984): Angenommen $a^p + b^p = c^p$ für Primzahl $p \geq 5$. Dann liefert die **Frey-Kurve** $y^2 = x(x - a^p)(x + b^p)$ eine semistabile elliptische Kurve.

2. **Ribet** (1990): Beweist Epsilon-Vermutung: Die Frey-Kurve ist **nicht modular** (falls sie existiert).

3. **Wiles** (1995): Beweist Shimura-Taniyama für alle semistabilen Kurven.

4. **Widerspruch**: Frey-Kurve müsste modular sein (Wiles) und gleichzeitig nicht modular (Ribet). Also kann $a^p + b^p = c^p$ keine Lösung haben.

---

## 8. Elliptic Curve Cryptography

### ECDH-Protokoll

Das **Elliptic Curve Diffie-Hellman** Protokoll:

**Öffentliche Parameter:** Kurve $E/\mathbb{F}_p$, Basispunkt $G$ der Ordnung $n$

**Schlüsselerzeugung:**
- Alice: Privat $a \in \{1, \ldots, n-1\}$, Öffentlich $A = aG$
- Bob: Privat $b \in \{1, \ldots, n-1\}$, Öffentlich $B = bG$

**Gemeinsames Geheimnis:**
$$S = aB = a(bG) = ab \cdot G = b(aG) = bA$$

**Sicherheit:** Basiert auf dem **Diskreten Logarithmusproblem** in $E(\mathbb{F}_p)$:

> Gegeben $P$ und $Q = kP$, finde $k$.

Bestes bekanntes Verfahren: Baby-Step Giant-Step in $\mathcal{O}(\sqrt{n})$ Zeit.

### Vergleich mit RSA

| Sicherheitsniveau | RSA Schlüssellänge | ECC Schlüssellänge |
|------------------|-------------------|-------------------|
| 80 Bit | 1024 Bit | 160 Bit |
| 128 Bit | 3072 Bit | 256 Bit |
| 256 Bit | 15360 Bit | 512 Bit |

ECC bietet gleiche Sicherheit bei deutlich kürzeren Schlüsseln.

---

## 9. ECM-Faktorisierung

### Lenstra's Elliptic Curve Method

**Idee:** Verwende die Gruppenstruktur von $E(\mathbb{F}_p)$ für verschiedene Primteiler $p$ von $n$.

**Algorithmus:**
1. Wähle zufällige Kurve $E$ und Punkt $P$ über $\mathbb{Z}/n\mathbb{Z}$
2. Berechne $k \cdot P$ für $k = \text{kgV}(1, 2, \ldots, B)$
3. Wenn $\gcd(\text{Nenner}, n) \neq 1$: Faktor gefunden!
4. Andernfalls: Neue Kurve wählen

**Warum funktioniert das?**

Sei $p | n$ ein Primteiler. In $E(\mathbb{F}_p)$ gilt:

$$\#E(\mathbb{F}_p) \approx p$$

Wenn $\#E(\mathbb{F}_p)$ **$B$-glatt** ist (alle Primfaktoren $\leq B$), dann ist $k$ ein Vielfaches von $\#E(\mathbb{F}_p)$, und damit $kP \equiv \mathcal{O} \pmod{p}$ aber nicht $\pmod{n}$.

**Laufzeit:** $L_p[1/2, \sqrt{2}]$ für den kleinsten Primteiler $p$ — optimal für mittlere Faktoren.

---

## 10. Kongruente Zahlen

### Definition

Eine positive ganze Zahl $n$ heißt **kongruente Zahl**, wenn sie der Flächeninhalt eines rechtwinkligen Dreiecks mit rationalen Seiten ist.

**Beispiele:**
- $n = 5$: kongruent (Seiten $\frac{3}{2}, \frac{20}{3}, \frac{41}{6}$)
- $n = 6$: kongruent (3-4-5-Dreieck, Fläche = 6)
- $n = 1$: nicht kongruent (Tunnell-Theorem, unter BSD bewiesen)

### Verbindung zur Elliptischen Kurve

Die **kongruente Zahl Kurve** ist:

$$E_n: y^2 = x^3 - n^2 x$$

**Satz** (Tunnell 1983, unter BSD): $n$ ist eine kongruente Zahl genau dann wenn $\text{Rang}(E_n(\mathbb{Q})) > 0$.

Dies verbindet ein 2000 Jahre altes geometrisches Problem mit der modernen Theorie elliptischer Kurven!

---

## 11. Bekannte Kurven

### secp256k1 (Bitcoin)

$$y^2 = x^3 + 7 \pmod{p}$$

$$p = 2^{256} - 2^{32} - 977$$

Eigenschaften:
- $a = 0$: Besonders effizient (kein $a$-Term in Formeln)
- Gruppenordnung $n$ ist eine Primzahl (~$2^{256}$)
- Basispunkt $G$ ist Teil des Standards SEC2

Verwendet in: Bitcoin, Ethereum (Signaturverfahren ECDSA)

### Curve25519

$$y^2 = x^3 + 486662 x^2 + x \pmod{p}$$

$$p = 2^{255} - 19$$

(Montgomery-Form, optimal für Diffie-Hellman)

Eigenschaften:
- Cofaktor 8 (nicht primzahlge Ordnung, dafür effizienter)
- Seitenkanalresistent durch konstante Laufzeit
- Gruppenordnung $= 8 \cdot$ Primzahl

Verwendet in: Signal-Protokoll, WireGuard, TLS 1.3, SSH

### y² = x³ − x (BSD-Beispiel)

- $a = -1$, $b = 0$
- $j = 1728$ (CM durch $\mathbb{Z}[i]$)
- $\text{Rang} = 0$ über $\mathbb{Q}$
- $E(\mathbb{Q})_{\text{tors}} = \mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/2\mathbb{Z}$
- Konsistent mit BSD: $L(E, 1) \neq 0$

---

## Supersingularität

Eine Kurve $E/\mathbb{F}_p$ heißt **supersingulär** wenn einer der äquivalenten Bedingungen gilt:

1. $\#E(\mathbb{F}_p) \equiv 1 \pmod{p}$
2. Die Frobenius-Spur $a_p \equiv 0 \pmod{p}$
3. Der $p$-Rang von $E$ ist 0
4. Der Endomorphismus-Ring ist nichtkommutativ (Quaternionenalgebra)

Supersingulare Kurven spielen eine Schlüsselrolle in der isogenie-basierten Kryptographie (SIDH, CSIDH) als Post-Quantum-Kandidaten.

---

*Dokumentation erstellt mit KaTeX-Unterstützung für mathematische Formeln.*
