# mersenne_fermat.py — Mersenne- und Fermat-Primzahlen

**Autor**: Michael Fuhrmann
**Letzte Änderung**: 2026-03-12
**Modul**: `src/mersenne_fermat.py`

---

## Übersicht

Dieses Modul implementiert zwei klassische Primzahlklassen der Zahlentheorie:
**Mersenne-Primzahlen** ($M_p = 2^p - 1$) und **Fermat-Primzahlen** ($F_n = 2^{2^n} + 1$).

---

## Klasse `MersennePrimes`

### Mathematischer Hintergrund

Eine **Mersenne-Zahl** $M_p = 2^p - 1$ ist eine Zahl der Form Zweierpotenz minus eins.
Notwendige Bedingung für Primalität: der Exponent $p$ muss selbst prim sein.

#### Lucas-Lehmer-Test (THEOREM, vollständig bewiesen)

Der effizienteste deterministische Test für Mersenne-Zahlen:

$$S_0 = 4, \quad S_k = S_{k-1}^2 - 2 \pmod{M_p}$$

**Satz**: Für primes $p \geq 3$ gilt: $M_p$ ist prim $\iff S_{p-2} \equiv 0 \pmod{M_p}$.

Laufzeit: $O(p^2 \log p)$ Bitoperationen (mit FFT-Multiplikation sogar $O(p \log^2 p)$).

#### Wagstaff-Heuristik (CONJECTURE)

Die Wahrscheinlichkeit, dass $M_p$ prim ist, beträgt näherungsweise:

$$\Pr[M_p \text{ prim}] \approx \frac{e^\gamma \cdot \ln 2}{\ln p}$$

wobei $\gamma \approx 0.5772...$ die **Euler-Mascheroni-Konstante** ist.

**Wichtig**: Dies ist eine unbewiese Heuristik (Wagstaff 1983). Sie impliziert (als Conjecture),
dass unendlich viele Mersenne-Primzahlen existieren.

#### Euklid-Euler-Theorem (THEOREM, vollständig bewiesen)

Eine gerade Zahl $n$ ist vollkommen genau dann, wenn:
$$n = 2^{p-1} \cdot (2^p - 1) \quad \text{mit } M_p \text{ prim}$$

- Euklid (~300 v. Chr.): Hinrichtung (wenn $M_p$ prim, dann $n$ vollkommen)
- Euler (1747): Rückrichtung (jede gerade vollkommene Zahl hat diese Form)

Ob **ungerade** vollkommene Zahlen existieren, ist ein **offenes Problem**.

### Bekannte Mersenne-Primzahlen (51 Stück, Stand 2024)

| Exponent $p$ | Stellenzahl von $M_p$ | Entdeckungsjahr |
|---|---|---|
| 2, 3, 5, 7 | 1, 1, 2, 3 | Antike |
| 13–31 | 4–10 | 1461–1772 |
| 61–521 | 19–157 | 1883–1952 |
| ... | ... | ... |
| 82589933 | ~24.8 Mio. | 2018 (GIMPS) |

### Methoden

| Methode | Beschreibung | Komplexität |
|---|---|---|
| `lucas_lehmer_test(p)` | Deterministischer Primzahltest für $M_p$ | $O(p^2)$ |
| `mersenne_number(p)` | Berechnet $2^p - 1$ | $O(p)$ |
| `wagstaff_heuristic(p)` | Heuristische Primwahrscheinlichkeit | $O(1)$ |
| `perfect_number_from_mersenne(p)` | Vollkommene Zahl $2^{p-1} \cdot M_p$ | $O(p^2)$ |
| `verify_perfect_number(n)` | Prüft $\sigma(n) = 2n$ | $O(\sqrt{n})$ |
| `get_first_n_mersenne_primes(n)` | Erste $n$ bekannte Exponenten | $O(1)$ |
| `digit_count(p)` | Anzahl Stellen von $M_p$ | $O(1)$ |

---

## Klasse `FermatPrimes`

### Mathematischer Hintergrund

**Fermat-Zahlen** haben die Form $F_n = 2^{2^n} + 1$:

| $n$ | $F_n$ | Status |
|---|---|---|
| 0 | 3 | **PRIM** |
| 1 | 5 | **PRIM** |
| 2 | 17 | **PRIM** |
| 3 | 257 | **PRIM** |
| 4 | 65537 | **PRIM** |
| 5 | 4294967297 | Zusammengesetzt: $641 \times 6700417$ (Euler 1732) |
| 6 | ~ $1.8 \times 10^{19}$ | Zusammengesetzt: $274177 \times 67280421310721$ |
| $\geq 5$ | riesig | Alle bekannten zusammengesetzt |

#### Pépin-Test (THEOREM, Pépin 1877)

Für $n \geq 1$ gilt: $F_n$ ist prim $\iff$:
$$3^{(F_n - 1)/2} \equiv -1 \pmod{F_n}$$

Basiert auf dem **Euler-Kriterium**: Für Primzahl $p$ und $\gcd(a,p)=1$ gilt
$a^{(p-1)/2} \equiv \left(\frac{a}{p}\right) \pmod{p}$ (Legendre-Symbol).

Da 3 ein quadratischer Nichtrest modulo $F_n$ ist (für primes $F_n$, $n \geq 1$),
ist $\left(\frac{3}{F_n}\right) = -1$.

#### Fermat-Heuristik (CONJECTURE)

Die Wahrscheinlichkeit, dass $F_n$ prim ist, beträgt:
$$\Pr[F_n \text{ prim}] \approx \frac{1}{2^n \cdot \ln 2}$$

Die Summe $\sum_{n=0}^{\infty} \Pr[F_n \text{ prim}]$ **konvergiert** —
was suggeriert, dass nur endlich viele Fermat-Primzahlen existieren.

**CONJECTURE**: Nur $F_0, \ldots, F_4$ sind prim. Nicht bewiesen.

#### Gauß-Wantzel-Theorem (BEWIESEN)

Ein reguläres $p$-Eck ist mit Zirkel und Lineal konstruierbar genau dann, wenn
$p$ ein Produkt aus einer Zweierpotenz und **verschiedenen Fermat-Primzahlen** ist.

Insbesondere: reguläres $F_n$-Eck konstruierbar $\iff$ $F_n$ prim.

### Methoden

| Methode | Beschreibung |
|---|---|
| `fermat_number(n)` | Berechnet $F_n = 2^{2^n} + 1$ |
| `pepin_test(n)` | Deterministischer Primalitätstest für $F_n$ |
| `is_prime_fermat(n)` | Kombinierter Test (bekannte Ergebnisse + Pépin) |
| `get_known_factor(n)` | Bekannter Faktor (falls vorhanden) |
| `verify_factorization(n)` | Verifiziert bekannte Faktorisierung |
| `wagstaff_heuristic_fermat()` | Heuristische Wahrscheinlichkeiten für $n=0..20$ |
| `gauss_wantzel_constructible(n)` | Konstruierbarkeit des $F_n$-Ecks |

---

## Beispielverwendung

```python
from mersenne_fermat import MersennePrimes, FermatPrimes

# Mersenne
mp = MersennePrimes()
print(mp.lucas_lehmer_test(61))        # True — M_61 ist prim
print(mp.perfect_number_from_mersenne(5))  # 496
print(mp.wagstaff_heuristic(100))      # ca. 0.12 (Heuristik!)

# Fermat
fp = FermatPrimes()
print(fp.pepin_test(4))                # True — F_4 = 65537 prim
print(fp.pepin_test(5))                # False — F_5 zusammengesetzt
print(fp.get_known_factor(5))          # 641
print(fp.verify_factorization(5))      # True
```

---

## Offene Probleme (Stand 2024)

1. **Unendlich viele Mersenne-Primzahlen?** (CONJECTURE — nicht bewiesen)
2. **Nur F_0–F_4 Fermat-Primzahlen?** (CONJECTURE — nicht bewiesen)
3. **Ungerade vollkommene Zahlen?** (Offen — keine gefunden, kein Unmöglichkeitsbeweis)
4. **GIMPS-Vollständigkeit**: Ist M_82589933 die 51. und letzte? (unbekannt)

---

## Quellen

- Lucas, É. (1878). *Théorie des fonctions numériques simplement périodiques*
- Lehmer, D.H. (1930). *An extended theory of Lucas functions*
- Pépin, T. (1877). *Sur la formule $2^{2^n} + 1$*
- Wagstaff, S.S. (1983). *Divisors of Mersenne numbers*
- GIMPS (2018). *Discovery of M_82589933*, www.mersenne.org
