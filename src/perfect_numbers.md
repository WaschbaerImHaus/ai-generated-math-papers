# perfect_numbers.py – Vollkommene Zahlen

**Autor:** Michael Fuhrmann
**Erstellt:** 2026-03-12
**Letzte Änderung:** 2026-03-12
**Build:** 122

---

## Überblick

Dieses Modul implementiert die Theorie der vollkommenen Zahlen (perfekten Zahlen),
eines der ältesten Themen der Zahlentheorie, das bis in die Antike zurückreicht.

---

## Definition

Eine natürliche Zahl $n$ heißt **vollkommen** (oder **perfekt**), wenn sie gleich
der Summe ihrer **echten** Teiler ist:

$$\sigma(n) = 2n$$

äquivalent: die Summe aller echten Teiler (ohne $n$ selbst) ist gleich $n$.

**Beispiele:**
- $6 = 1 + 2 + 3$ ✓
- $28 = 1 + 2 + 4 + 7 + 14$ ✓
- $12 \neq 1+2+3+4+6 = 16$ ✗ (abundant)

---

## σ-Funktion

Die **Teiler-Summenfunktion** $\sigma(n) = \sum_{d \mid n} d$ ist zentral für die
Theorie vollkommener Zahlen.

### Eigenschaften

**Multiplikativität:** Für $\gcd(m, n) = 1$:
$$\sigma(mn) = \sigma(m) \cdot \sigma(n)$$

**Primzahlpotenzen:**
$$\sigma(p^k) = 1 + p + p^2 + \cdots + p^k = \frac{p^{k+1} - 1}{p - 1}$$

### Klassifikation

| Typ | Bedingung | Abundanzindex $\sigma(n)/n$ |
|-----|-----------|---------------------------|
| Vollkommen | $\sigma(n) = 2n$ | $= 2$ |
| Abundant | $\sigma(n) > 2n$ | $> 2$ |
| Defizient | $\sigma(n) < 2n$ | $< 2$ |

---

## Gerade vollkommene Zahlen

### Euklid-Euler-Theorem (vollständig bewiesen)

$$n \text{ gerade vollkommen} \iff n = 2^{p-1}(2^p - 1) \text{ mit } 2^p - 1 \text{ prim}$$

- **Euklid (~300 v. Chr.):** Bewies die Hinrichtung (Suffizienzbedingung)
- **Euler (1849):** Bewies die Rückrichtung (Notwendigkeitsbedingung)

Zahlen der Form $M_p = 2^p - 1$ mit $p$ prim heißen **Mersenne-Zahlen**.
Ist $M_p$ zusätzlich prim, heißt sie **Mersenne-Primzahl**.

### Erste gerade vollkommene Zahlen

| $p$ | $M_p = 2^p - 1$ | $n = 2^{p-1} M_p$ | Stellen |
|-----|-----------------|-------------------|---------|
| 2   | 3               | **6**             | 1       |
| 3   | 7               | **28**            | 2       |
| 5   | 31              | **496**           | 3       |
| 7   | 127             | **8128**          | 4       |
| 13  | 8191            | **33550336**      | 8       |
| 17  | 131071          | **8589869056**    | 10      |
| 19  | 524287          | **137438691328**  | 12      |
| 31  | 2147483647      | $2^{30}(2^{31}-1)$| 19      |

### Offene Fragen (CONJECTURE)

- Gibt es **unendlich viele** gerade vollkommene Zahlen?
  (Äquivalent: gibt es unendlich viele Mersenne-Primzahlen?)
- Stand 2024: **51 bekannte** gerade vollkommene Zahlen

---

## Ungerade vollkommene Zahlen

### Hauptvermutung

> **CONJECTURE (seit Antike, offen):** Es gibt keine ungeraden vollkommenen Zahlen.

Trotz aller Bemühungen ist weder ein Beweis noch ein Gegenbeispiel bekannt.

### Bekannte notwendige Bedingungen (bewiesen)

Falls $n$ eine ungerade vollkommene Zahl ist, dann gilt zwingend:

**1. Euler-Form (Euler ~1849):**
$$n = p^{4a+1} \cdot m^2, \quad p \equiv 1 \pmod{4}, \quad p \text{ prim}, \quad \gcd(p, m) = 1$$

**2. Größenschranke:**
$$n > 10^{1500} \quad \text{[aktuelle Schranke, 2020er]}$$
(Brent-Cohen-te Riele 1991 zeigten $n > 10^{300}$, seitdem verschärft)

**3. Nielsen-Schranke (2006):**
$$\omega(n) \geq 9 \quad \text{(mindestens 9 verschiedene Primfaktoren)}$$

**4. Chein-Schranke (1979):**
$$\Omega(n) \geq 101 \quad \text{(mindestens 101 Primfaktoren mit Vielfachheit)}$$

**5. Touchard-Kongruenz (1953):**
$$n \equiv 1 \pmod{12} \quad \text{ODER} \quad n \equiv 9 \pmod{36}$$

### Kombination der Bedingungen

Eine Zahl, die eine dieser notwendigen Bedingungen verletzt, **kann nicht**
ungerade und vollkommen sein. Alle bekannten kleinen ungeraden Zahlen
werden durch den σ-Test oder die Touchard-Kongruenz ausgeschlossen.

---

## Klassen

### `SigmaFunction`

| Methode | Beschreibung |
|---------|-------------|
| `sigma(n, k)` | Berechnet $\sigma_k(n)$ |
| `sigma_from_factorization(dict)` | $\sigma$ aus Primfaktorzerlegung |
| `is_perfect(n)` | Prüft $\sigma(n) = 2n$ |
| `is_abundant(n)` | Prüft $\sigma(n) > 2n$ |
| `is_deficient(n)` | Prüft $\sigma(n) < 2n$ |
| `abundancy_index(n)` | Berechnet $\sigma(n)/n$ |
| `multiplicativity_check(m, n)` | Verifiziert $\sigma(mn) = \sigma(m)\sigma(n)$ |
| `perfect_numbers_up_to(limit)` | Findet alle vollkommenen Zahlen $\leq$ limit |

### `EvenPerfectNumbers`

| Methode | Beschreibung |
|---------|-------------|
| `mersenne_number(p)` | Berechnet $M_p = 2^p - 1$ |
| `is_mersenne_prime(p)` | Prüft ob $M_p$ prim ist |
| `even_perfect_from_exponent(p)` | Berechnet $2^{p-1}(2^p-1)$ |
| `first_n_even_perfect(n)` | Erste $n$ gerade vollkommene Zahlen |
| `verify_euclid_euler(p)` | Vollständige Verifikation des Theorems |
| `even_perfect_digit_count(p)` | Anzahl Dezimalstellen |

### `OddPerfectNumberBounds`

| Methode | Beschreibung |
|---------|-------------|
| `euler_form_check(n)` | Prüft Euler-Form $n = p^{4a+1}m^2$ |
| `touchard_congruence(n)` | Prüft Touchard $n \equiv 1 \pmod{12}$ oder $\equiv 9 \pmod{36}$ |
| `nielsen_bound_check(n)` | Prüft Nielsen $\omega(n) \geq 9$ |
| `chein_structure(n)` | Analysiert Chein-Exponenten-Struktur |
| `is_excluded_as_odd_perfect(n)` | Kombinierter Ausschlusstest |
| `summary_necessary_conditions()` | Zusammenfassung aller Bedingungen |

---

## Verwendungsbeispiele

```python
from perfect_numbers import SigmaFunction, EvenPerfectNumbers, OddPerfectNumberBounds

# Sigma-Funktion
print(SigmaFunction.sigma(6))    # 12 = 1+2+3+6
print(SigmaFunction.is_perfect(6))  # True (σ(6)=12=2·6)
print(SigmaFunction.abundancy_index(12))  # > 2 (abundant)

# Gerade vollkommene Zahlen
epn = EvenPerfectNumbers()
print(epn.first_n_even_perfect(5))
# [6, 28, 496, 8128, 33550336]

# Verifikation Euklid-Euler für p=5
result = EvenPerfectNumbers.verify_euclid_euler(5)
print(result)
# {'exponent_p': 5, 'p_is_prime': True, 'mersenne_number': 31,
#  'mersenne_is_prime': True, 'perfect_number': 496, 'sigma_equals_2n': True}

# Touchard-Kongruenz
tc = OddPerfectNumberBounds.touchard_congruence(6)
print(tc['satisfies_touchard'])  # False → 6 kann nicht ungerade vollkommen sein

# Ausschlusstest
excluded, reasons = OddPerfectNumberBounds.is_excluded_as_odd_perfect(15)
print(excluded, reasons)  # True, [...]
```

---

## Historische Anmerkungen

- **Euklid (~300 v. Chr.):** Erste 4 vollkommene Zahlen bekannt (6, 28, 496, 8128)
- **Nicomachus (100 n. Chr.):** Erste 4 vollkommene Zahlen katalogisiert
- **Euler (1707–1783):** Vollständige Charakterisierung gerader vollkommener Zahlen
- **Touchard (1953):** Kongruenzschranken für ungerade vollkommene Zahlen
- **Brent, Cohen, te Riele (1991):** Untere Schranke $n > 10^{300}$

---

## Offene Probleme

1. **Existenz ungerader vollkommener Zahlen** (CONJECTURE, seit Antike, offen)
2. **Unendlichkeit der geraden vollkommenen Zahlen** (äquivalent zu unendlich vielen Mersenne-Primzahlen, CONJECTURE)
3. **Verbesserung der Nielsen-Schranke** (kann man $\omega(n) \geq 10$ oder mehr zeigen?)

---

## Literatur

- Euklid, *Elemente*, Buch IX, Satz 36 (~300 v. Chr.)
- Euler, L. (1849). *Opera Postuma.* (Charakterisierung gerader vollkommener Zahlen)
- Touchard, J. (1953). On prime numbers and perfect numbers. *Scripta Math.* 19.
- Nielsen, P.P. (2006). Odd perfect numbers have at least nine distinct prime factors. *Math. Comp.* 76.
- Brent, R.P., Cohen, G.L., te Riele, H.J.J. (1991). Improved techniques for lower bounds for odd perfect numbers. *Math. Comp.* 57.
