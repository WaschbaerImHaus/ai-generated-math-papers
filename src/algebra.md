# algebra.py – Algebra-Modul

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/algebra.py`

---

## Überblick

Dieses Modul implementiert grundlegende und fortgeschrittene algebraische Konzepte:
- Polynomklasse mit vollständigen Operationen
- Gleichungslöser für lineare und quadratische Gleichungen
- Zahlentheoretische Funktionen (ggT, kgV, Primzahlen, Euler-Phi)

---

## Klasse: `Polynomial`

Repräsentiert ein Polynom der Form:

```
p(x) = a_n·xⁿ + a_{n-1}·x^{n-1} + ... + a₁·x + a₀
```

Die Koeffizienten werden **von höchstem Grad zu niedrigstem** gespeichert:

```python
p = Polynomial([1, -3, 2])   # repräsentiert: x² - 3x + 2
```

### Konstruktor

```python
Polynomial(coefficients: list)
```

Führende Nullen werden automatisch entfernt (Normalisierung).

### Eigenschaften

| Eigenschaft | Beschreibung |
|-------------|-------------|
| `degree` | Grad des Polynoms (höchster Exponent) |
| `coefficients` | Liste der Koeffizienten (intern) |

### Methoden

#### `evaluate(x: float) → float`
Wertet das Polynom an der Stelle `x` aus.
**Algorithmus:** Horner-Schema – reduziert Multiplikationen von O(n²) auf O(n):
```
a·x³ + b·x² + c·x + d = ((a·x + b)·x + c)·x + d
```

#### `derivative() → Polynomial`
Berechnet die symbolische Ableitung via Potenzregel:
```
d/dx (a·xⁿ) = a·n·xⁿ⁻¹
```

#### `__add__`, `__sub__`, `__mul__`
Addition, Subtraktion und Multiplikation von Polynomen.
Multiplikation entspricht einer **Faltungsoperation**:
```
c_k = Σ aᵢ·bⱼ  für alle i+j=k
```

#### `__str__() → str`
Gibt das Polynom in mathematischer Notation aus (z.B. `"x^2 - 3x + 2"`).

---

## Funktionen: Gleichungslöser

### `solve_linear(a, b) → float`
Löst `ax + b = 0` → `x = -b/a`.

**Sonderfälle:**
- `a = 0, b ≠ 0` → `ValueError` (kein Lösung)
- `a = 0, b = 0` → `ValueError` (unendlich viele Lösungen)

### `solve_quadratic(a, b, c) → list`
Löst `ax² + bx + c = 0` via **Mitternachtsformel (pq-Formel)**:

```
x₁,₂ = (-b ± √(b² - 4ac)) / (2a)
```

**Diskriminante D = b² - 4ac:**
| D | Ergebnis |
|---|----------|
| D > 0 | Zwei verschiedene reelle Wurzeln |
| D = 0 | Eine doppelte reelle Wurzel |
| D < 0 | Zwei komplexkonjugierte Wurzeln |

---

## Funktionen: Zahlentheorie

### `gcd(a, b) → int`
Größter gemeinsamer Teiler (ggT) via **Euklidischen Algorithmus**:

```
ggT(a, b) = ggT(b, a mod b)
```

**Laufzeit:** O(log min(a,b))

**Beispiel:**
```
ggT(48, 18) = ggT(18, 12) = ggT(12, 6) = ggT(6, 0) = 6
```

### `lcm(a, b) → int`
Kleinstes gemeinsames Vielfaches (kgV) via:
```
kgV(a, b) = |a| // ggT(a,b) × |b|
```
(Frühzeitige Division verhindert Integer-Überlauf)

### `extended_gcd(a, b) → tuple`
Erweiterter Euklidischer Algorithmus – findet `x, y` mit:
```
a·x + b·y = ggT(a, b)   (Bezout-Gleichung)
```

**Anwendung:** Modulare Inverse, RSA-Kryptographie

### `mod_inverse(a, m) → int`
Berechnet `x` mit `(a·x) % m == 1`.
Existiert genau dann wenn `ggT(a, m) = 1`.

### `is_prime(n) → bool`
**Primzahltest durch Probedivision** mit drei Optimierungen:
1. Sonderfälle n < 4
2. Teiler 2 und 3 ausschließen
3. Nur Zahlen der Form `6k ± 1` testen bis √n

**Laufzeit:** O(√n)

### `prime_factorization(n) → dict`
Primfaktorzerlegung: `n = p₁^e₁ · p₂^e₂ · ... · pₖ^eₖ`

**Rückgabe:** `{primzahl: exponent}`, z.B. `{2: 3, 3: 2, 5: 1}` für 360

### `euler_phi(n) → int`
**Eulers Phi-Funktion:** Anzahl der Zahlen von 1 bis n, die zu n teilerfremd sind.

**Berechnung via Primfaktorzerlegung:**
```
φ(n) = n · ∏(1 - 1/p)   für alle Primfaktoren p von n
```

**Wichtige Eigenschaften:**
- φ(1) = 1
- φ(p) = p-1 für Primzahlen p
- φ(pᵏ) = p^{k-1}·(p-1)
- φ ist multiplikativ: φ(m·n) = φ(m)·φ(n) wenn ggT(m,n)=1

**Anwendung:** Eulers Satz (`a^φ(n) ≡ 1 mod n`), RSA-Verschlüsselung

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | Quadratwurzeln, sqrt |
| `cmath` | Komplexe Quadratwurzeln (für quadratische Gleichungen) |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from algebra import Polynomial, solve_quadratic, gcd, euler_phi

# Polynom anlegen und auswerten
p = Polynomial([1, -3, 2])   # x² - 3x + 2
print(p)               # "x^2 - 3x + 2"
print(p.evaluate(1))   # 0.0  (Nullstelle)
print(p.derivative())  # "2x - 3"

# Quadratische Gleichung lösen
roots = solve_quadratic(1, -3, 2)
print(roots)   # [2.0, 1.0]

# Zahlentheorie
print(gcd(48, 18))       # 6
print(euler_phi(12))     # 4  (4 Zahlen ≤12 teilerfremd zu 12: 1,5,7,11)
```

---

## Mathematischer Hintergrund

### Euklidischer Algorithmus (ca. 300 v. Chr.)
Einer der ältesten Algorithmen der Mathematik. Basiert auf der Beobachtung, dass `ggT(a,b) = ggT(b, a mod b)`, da jeder gemeinsame Teiler von `a` und `b` auch ein Teiler von `a mod b` ist.

### Bezout-Identität
Für ganze Zahlen `a` und `b` existieren immer ganze Zahlen `x, y` mit `ax + by = ggT(a,b)`. Der erweiterte Euklidische Algorithmus findet diese Koeffizienten in O(log min(a,b)) Schritten.

### Horner-Schema (William George Horner, 1819)
Effiziente Polynomauswertung durch verschachtelte Multiplikation. Vermeidet wiederholte Potenzberechnungen und benötigt nur n Multiplikationen und n Additionen für ein Polynom vom Grad n.
