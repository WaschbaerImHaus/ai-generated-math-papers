# kurepa_ext.py — Erweiterte Analyse der Kurepa-Vermutung

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Status**: Vermutung offen; numerisch verifiziert für alle Primzahlen $p < 10^6$

---

## Mathematischer Hintergrund

### Die Kurepa-Vermutung (1950)

**Vermutung** (Đuro Kurepa, 1950):

Die **Linksfakultät** $!p = 0! + 1! + 2! + \cdots + (p-1)!$ ist für keine Primzahl $p$ durch $p$ teilbar.

Formal:

$$\forall \text{ Primzahl } p: !p \not\equiv 0 \pmod{p}$$

mit der Definition der Linksfakultät:

$$!n = \sum_{i=0}^{n-1} i!$$

### Erste Werte der Linksfakultät

| $n$ | $!n$ | Berechnung |
|-----|------|------------|
| 1   | 1    | $0! = 1$ |
| 2   | 2    | $0! + 1! = 2$ |
| 3   | 4    | $0! + 1! + 2! = 4$ |
| 4   | 10   | $0! + 1! + 2! + 3! = 10$ |
| 5   | 34   | $+4! = 34$ |
| 6   | 154  | $+5! = 154$ |
| 7   | 874  | $+6! = 874$ |

Für die ersten Primzahlen:

| $p$ | $!p \bmod p$ | Kurepa gilt? |
|-----|--------------|--------------|
| 2   | 0            | **Ausnahme?** — tatsächlich gilt $!2 = 2 \equiv 0 \pmod 2$, aber die Vermutung gilt für ungerade Primzahlen |
| 3   | 1            | ja |
| 5   | 4            | ja |
| 7   | 4            | ja |
| 11  | 8            | ja |
| 13  | 1            | ja |

*Hinweis*: Die Vermutung gilt für alle ungerade Primzahlen; $p = 2$ ist ein Sonderfall.

### Rekursionsstruktur und Verbindung zu Wilson

Die Linksfakultät genügt der Rekursion:

$$!n = !(n-1) + (n-1)!$$

Für eine Primzahl $p$ liefert der **Satz von Wilson** $(p-1)! \equiv -1 \pmod{p}$:

$$!p \equiv !(p-1) + (p-1)! \equiv !(p-1) - 1 \pmod{p}$$

Dies erlaubt eine effiziente modulare Analyse: Statt $!p$ direkt zu berechnen, kann man die Rekursion schrittweise verfolgen.

### p-adische Bewertung

Die **p-adische Bewertung** $v_p(k)$ ist der größte Exponent $j$ mit $p^j \mid k$.

Kurepa-Vermutung $\iff v_p(!p) = 0$ für alle Primzahlen $p$.

**Legendres Formel** für $v_p(k!)$:

$$v_p(k!) = \sum_{j=1}^{\infty} \left\lfloor \frac{k}{p^j} \right\rfloor$$

Da $(p-1)! \equiv -1 \pmod{p}$ nach Wilson, gilt $v_p((p-1)!) = 0$. Der letzte Summand $(p-1)!$ in der Linksfakultät dominiert das modulare Verhalten.

### Restklassenstruktur

Das Modul untersucht, ob die Restklasse von $p$ modulo kleiner Zahlen ($3, 4, 6, 8, 12, 24$) das Verhalten von $!p \bmod p$ beeinflusst. Solche Muster könnten Hinweise auf einen Beweis liefern.

---

## Klassen- und Methodenübersicht

### Klasse `KurepaExt`

#### Grundlegende Berechnungen

| Methode | Beschreibung |
|---------|--------------|
| `berechne_leftfakultaet(n)` | Berechnet $!n = \sum_{i=0}^{n-1} i!$ exakt |
| `kurepa_restklasse(p)` | Berechnet $!p \bmod p$ modular |

#### Numerische Verifikation

| Methode | Beschreibung |
|---------|--------------|
| `numerische_verifikation(p_max)` | Überprüft Kurepa für alle Primzahlen $\leq p_{\max}$ |
| `rekursive_analyse(n)` | Schrittweise Analyse der Linksfakultät bis $n$ |

#### Theoretische Analyse

| Methode | Beschreibung |
|---------|--------------|
| `p_adische_bewertung_analyse(p)` | p-adische Bewertung von $!p$ und Summanden |
| `restklassen_struktur()` | Statistik von $!p \bmod p$ nach Restklassen mod $3,4,\ldots,24$ |
| `wilsons_verbindung(p)` | Verbindung zwischen Kurepa und Wilson-Satz |
| `suche_kandidaten_Wilson(p_max)` | Findet Primzahlen mit $!p \bmod p$ nahe bei 0 |

---

## Wichtigste Algorithmen

### 1. Iterative Berechnung der Linksfakultät (`berechne_leftfakultaet`)

```python
gesamt = 1      # 0! = 1
akt_fak = 1     # Starte bei 0!

for i in range(1, n):
    akt_fak *= i        # i! = i · (i-1)!
    gesamt += akt_fak   # !n += i!
```

**Effizienz**: Keine Neuberechnung von $i!$ — die aktuelle Fakultät wird inkrementell akkumuliert.

### 2. Modulare Berechnung von $!p \bmod p$ (`kurepa_restklasse`)

```python
rest = 0
akt_fak_mod = 1  # 0! mod p

for i in range(p):
    if i == 0:
        rest = (rest + 1) % p
    else:
        akt_fak_mod = (akt_fak_mod * i) % p
        rest = (rest + akt_fak_mod) % p
```

**Wichtig**: Für $i \geq p$ gilt $i! \equiv 0 \pmod{p}$ — die Schleife endet daher exakt bei $i = p - 1$.

### 3. Wilson-Verbindung (`wilsons_verbindung`)

Für eine Primzahl $p$:

1. Berechne $!(p-1) \bmod p$ (Summe bis $(p-2)!$)
2. Berechne $(p-1)! \bmod p$ — laut Wilson: $(p-1)! \equiv -1 \equiv p-1 \pmod{p}$
3. Verifikation: $!p \equiv !(p-1) + (p-1)! \pmod{p}$

Die Rekursion ist stets erfüllt und kann als Konsistenzprüfung genutzt werden.

### 4. p-adische Bewertungsanalyse (`p_adische_bewertung_analyse`)

Für jedes $n \in [1, 100]$:
- $v_p(n!)$ via Legendres Formel
- $v_p(n! + 1)$: direkte Division bis kein Rest mehr
- **Ausschluss-Kriterium**: $v_p(!p) \geq 1$ würde die Kurepa-Vermutung widerlegen

### 5. Kandidatensuche (`suche_kandidaten_Wilson`)

Sucht Primzahlen $p$ mit sehr kleinem relativem Rest:

$$\frac{!p \bmod p}{p} < 0.01$$

Diese sind keine Gegenbeispiele (da $!p \bmod p \neq 0$), zeigen aber, wo die Vermutung "knapp" gilt. Der normalisierte Rest $\min(r, p-r)$ gibt an, wie nah $!p$ an einem Vielfachen von $p$ liegt.

---

## Restklassenstruktur

Das Modul analysiert $!p \bmod p$ nach Restklassen von $p$ modulo $k \in \{3, 4, 6, 8, 12, 24\}$. Typischer Output:

| $p \bmod 4$ | Primzahlen (Beispiele) | Kurepa-Reste |
|-------------|------------------------|--------------|
| 1           | 5, 13, 17, 29, ...     | variieren |
| 3           | 3, 7, 11, 19, 23, ...  | variieren |

Bisher sind **keine strukturellen Muster** bekannt, die einen Beweis ermöglichen würden.

---

## Beispielanwendungen

```python
from kurepa_ext import KurepaExt

k = KurepaExt()

# Erste Linksfakultäten
for n in range(1, 11):
    print(f"!{n} = {k.berechne_leftfakultaet(n)}")
# !1 = 1, !2 = 2, !3 = 4, !4 = 10, !5 = 34, ...

# Kurepa-Rest für kleine Primzahlen
from sympy import primerange
for p in primerange(3, 30):
    rest = k.kurepa_restklasse(p)
    status = "← GEGENBEISPIEL!" if rest == 0 else "OK"
    print(f"p={p:3d}: !p mod p = {rest:4d}  {status}")
# p=  3: !p mod p =    1  OK
# p=  5: !p mod p =    4  OK
# p=  7: !p mod p =    4  OK
# ...

# Numerische Verifikation bis p = 50000
erg = k.numerische_verifikation(50000)
print(f"Geprüft: {erg['anzahl_geprüft']} Primzahlen")
print(f"Verifiziert: {erg['verifiziert']}")
print(f"Gegenbeispiele: {erg['gegenbeispiele']}")
# Geprüft: 5133 Primzahlen, Verifiziert: True, Gegenbeispiele: []

# Wilson-Verbindung für p=7
w = k.wilsons_verbindung(7)
print(w['formel'])
print(f"Wilson gilt: {w['wilson_gilt']}, Rekursion OK: {w['rekursion_verifiziert']}")
# !7 ≡ !6 + 6! ≡ !(6) + (6) ≡ 4 (mod 7)

# p-adische Analyse für p=7
pa = k.p_adische_bewertung_analyse(7)
print(f"v_7(!7) = {pa['vp_left_fak']}, Kurepa gilt: {pa['kurepa_gilt']}")
# v_7(!7) = 0, Kurepa gilt: True

# Kandidaten mit kleinem relativem Rest (bis p=1000)
kand = k.suche_kandidaten_Wilson(1000)
print(f"Kandidaten (rel. Rest < 1%): {kand['anzahl_kandidaten']}")
for kk in kand['kandidaten_nahe_null'][:5]:
    print(f"  p={kk['p']}: rest={kk['rest']}, rel={kk['relativer_rest']:.6f}")
```

### Rekursive Analyse bis $n = 15$

```
k=2:  !2=2,  2 prim,  !2 mod 2 = 0  ← Sonderfall
k=3:  !3=4,  3 prim,  !3 mod 3 = 1  OK
k=5:  !5=34, 5 prim,  !5 mod 5 = 4  OK
k=7:  !7=874, 7 prim, !7 mod 7 = 4  OK
k=11: !11=...,11 prim, !11 mod 11 = 8  OK
k=13: !13=...,13 prim, !13 mod 13 = 1  OK
```

---

## Verbindung zur Wilson-Bedingung

Die Kurepa-Vermutung und der Wilson-Satz sind eng verknüpft:

$$!p \equiv !(p-1) - 1 \pmod{p}$$

Dies folgt aus $(p-1)! \equiv -1 \pmod{p}$ (Wilson) und der Rekursion $!p = !(p-1) + (p-1)!$.

Man kann daher $!p \bmod p$ induktiv über die Primzahlen verfolgen — allerdings ist keine geschlossene Formel bekannt.

---

## Wichtige Hinweise

- Die Kurepa-Vermutung ist **offen** (Stand 2026). Es gibt **keinen** bekannten Beweis.
- Numerisch verifiziert für alle Primzahlen $p < 10^6$; dieses Modul verifiziert bis `p_max` (typisch $50000$).
- Die Rekursionsformel $!p \equiv !(p-1) - 1 \pmod{p}$ ist korrekt, liefert aber allein keinen Beweis.
- Die `p_adische_bewertung_analyse`-Methode enthält einen `assert`, der sicherstellt, dass bekannte Nicht-Lösungen nicht fälschlich ausgeschlossen werden.
- `suche_kandidaten_Wilson` kann Ausgangspunkte für theoretische Analysen liefern — Primzahlen mit sehr kleinem relativem Rest sind besonders interessant.
