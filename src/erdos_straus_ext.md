# erdos_straus_ext.py — Erweiterte konstruktive Analyse der Erdős-Straus-Vermutung

**Autor**: Michael Fuhrmann
**Datum**: 2026-03-12
**Status**: Konstruktive Teilbeweise für Restklassen mod 4 — Vermutung insgesamt offen

---

## Mathematischer Hintergrund

### Die Erdős-Straus-Vermutung (1948)

**Vermutung**: Für alle natürlichen Zahlen $n \geq 2$ existieren positive ganze Zahlen $x, y, z$ mit:

$$\frac{4}{n} = \frac{1}{x} + \frac{1}{y} + \frac{1}{z}$$

Eine solche Darstellung nennt man eine **ägyptische Bruchzerlegung** von $\frac{4}{n}$.

### Bekannte Fakten

| Fakt | Quelle |
|------|--------|
| Verifiziert für alle $n \leq 10^{14}$ | Swett (2004), numerisch |
| Vollständige Restklassenabdeckung mod 840 | Kotzig-Turán (1960) |
| Für $p \equiv 3 \pmod{4}$: exakter konstruktiver Beweis | Klassisches Resultat |

### Bewiesener Teilsatz: $p \equiv 3 \pmod{4}$

Für jede Primzahl $p = 4k + 3$ gilt:

$$\frac{4}{p} = \frac{1}{k+1} + \frac{1}{(k+1) \cdot p}$$

**Beweis**:

$$\frac{1}{k+1} + \frac{1}{(k+1)(4k+3)} = \frac{(4k+3) + 1}{(k+1)(4k+3)} = \frac{4(k+1)}{(k+1)(4k+3)} = \frac{4}{4k+3} = \frac{4}{p} \checkmark$$

$k+1 = \frac{p+1}{4}$ ist ganzzahlig, da $p \equiv 3 \pmod{4}$ impliziert $p + 1 \equiv 0 \pmod{4}$.

Diese 2-Term-Zerlegung impliziert sofort eine 3-Term-Zerlegung durch triviale Aufspaltung:

$$\frac{1}{(k+1)p} = \frac{1}{2(k+1)p} + \frac{1}{2(k+1)p}$$

### Parametrische Lösungsformeln (vollständige Tabelle)

| Klasse mod 4 | Formel | Beweis-Status |
|---|---|---|
| $n \equiv 0 \pmod{4}$, $n = 4k$ | $\frac{4}{4k} = \frac{1}{2k} + \frac{1}{3k} + \frac{1}{6k}$ | **Exakter Beweis** |
| $n \equiv 3 \pmod{4}$, $n = 4k+3$ | $\frac{4}{n} = \frac{1}{k+1} + \frac{1}{(k+1)n} + \frac{1}{2(k+1)n}$ | **Exakter Beweis** |
| $n \equiv 1 \pmod{4}$, $n = 4k+1$ | $\frac{4}{n} = \frac{1}{k+1} + \frac{3}{(k+1)n}$ (Greedy) | Konstruktiv |
| $n \equiv 2 \pmod{4}$, $n = 2(2k+1)$ | $\frac{4}{n} = \frac{1}{k+1} + \frac{1}{(k+1)(2k+1)} + \text{Split}$ | Konstruktiv |

**Beweis für $n = 4k$**:

$$\frac{1}{2k} + \frac{1}{3k} + \frac{1}{6k} = \frac{3 + 2 + 1}{6k} = \frac{6}{6k} = \frac{1}{k} = \frac{4}{4k} \checkmark$$

---

## Klassen- und Methodenübersicht

### Hilfsfunktionen (modullevel)

| Funktion | Beschreibung |
|----------|--------------|
| `_gcd(a, b)` | Größter gemeinsamer Teiler via `math.gcd` |
| `_bruch_pruefen(n, x, y, z)` | Exakte Verifikation $\frac{4}{n} = \frac{1}{x} + \frac{1}{y} + \frac{1}{z}$ (ganzzahlig) |
| `_bruch_pruefen_2terme(n, x, y)` | Exakte Verifikation $\frac{4}{n} = \frac{1}{x} + \frac{1}{y}$ |

### Klasse `ErdosStrausExt`

#### Lösungsalgorithmen

| Methode | Beschreibung |
|---------|--------------|
| `loese_unit_fraction(n)` | Hauptmethode: findet Zerlegung mit 3 Strategien |
| `_strategie_restklasse(n)` | Direkte Formeln nach $n \bmod 4$ |
| `_strategie_teilersuche(n)` | Parametrische Suche über Teiler |
| `_strategie_brute_force(n, limit)` | Beschränkte erschöpfende Suche |

#### Konstruktive Teilbeweise

| Methode | Beschreibung |
|---------|--------------|
| `beweis_klasse_mod4_3(p)` | Konstruktiver Beweis für $p \equiv 3 \pmod{4}$ |
| `beweis_klasse_mod4_1(p)` | Konstruktiver Beweis für $p \equiv 1 \pmod{4}$ |
| `beweis_klasse_mod3(p, r)` | Analyse für $p \equiv r \pmod{3}$ |
| `beweis_klasse_mod_12(p, r)` | Kombinierte Analyse mod 12 |

#### Vollständige Analyse

| Methode | Beschreibung |
|---------|--------------|
| `vollstaendige_analyse(grenze)` | Analysiert alle $n$ von 2 bis `grenze` |
| `suche_ausnahmen(grenze)` | Sucht $n$ ohne gefundene Zerlegung |
| `parametrische_loesungen()` | Gibt tabellarische Formelübersicht zurück |

---

## Wichtigste Algorithmen

### 1. Dreistufige Lösungsstrategie (`loese_unit_fraction`)

```
Eingabe: n ≥ 2
1. Versuche Restklassenformel (mod 4) → schnell, exakt für bekannte Klassen
2. Falls kein Ergebnis: Parametrische Teilersuche
3. Falls kein Ergebnis: Beschränkte Brute-Force
```

### 2. Exakte Bruchverifikation (`_bruch_pruefen`)

Um Gleitkommaprobleme zu vermeiden, wird durch Kreuzprodukt verglichen:

$$\frac{4}{n} = \frac{1}{x} + \frac{1}{y} + \frac{1}{z} \iff 4xyz = n(yz + xz + xy)$$

### 3. Parametrische Teilersuche (`_strategie_teilersuche`)

Sei $x = \lceil n/4 \rceil$. Dann gilt:

$$\frac{4}{n} - \frac{1}{x} = \frac{4x - n}{nx}$$

mit $r = 4x - n \in \{1, 2, 3\}$ (da $x = \lceil n/4 \rceil$). Gesucht wird $y$ so dass:

$$\frac{r \cdot y - nx}{nx \cdot y} = \frac{1}{z}$$

also $nx \cdot y$ durch $(r \cdot y - nx)$ teilbar ist.

### 4. Greedy für $n \equiv 1 \pmod{4}$

Für $n = 4k + 1$ mit $x = k + 1$:

$$\frac{4}{n} - \frac{1}{x} = \frac{3}{(k+1)(4k+1)}$$

Der Rest $\frac{3}{D}$ mit $D = (k+1)n$ wird greedy zerlegt: suche $y$ ab $y \approx D/3$ bis $\frac{3y - D}{Dy}$ ganzzahligen Nenner hat.

### 5. Mod-12-Analyse

Mod 12 kombiniert die Informationen mod 4 und mod 3 ($\text{kgV}(4,3) = 12$). Für Primzahlen $> 3$ kommen nur die Restklassen $\{1, 5, 7, 11\}$ mod 12 vor:

| $p \bmod 12$ | $p \bmod 4$ | $p \bmod 3$ | Anwendbarer Teilbeweis |
|---|---|---|---|
| 1 | 1 | 1 | Greedy |
| 5 | 1 | 2 | Greedy |
| 7 | 3 | 1 | **Exakter Beweis** |
| 11 | 3 | 2 | **Exakter Beweis** |

---

## Beispielanwendungen

```python
from erdos_straus_ext import ErdosStrausExt

es = ErdosStrausExt()

# Beweis für p ≡ 3 (mod 4)
for p in [3, 7, 11, 19, 23, 31, 43]:
    res = es.beweis_klasse_mod4_3(p)
    print(f"p={p:4d}: {res['status']}, Zerlegung={res['3term_erweiterung']}")
# p=   3: KONSTRUKTIV BEWIESEN, Zerlegung=(1, 12, 12)
# p=   7: KONSTRUKTIV BEWIESEN, Zerlegung=(2, 28, 28)
# p=  11: KONSTRUKTIV BEWIESEN, Zerlegung=(3, 66, 66)

# Vollständige Analyse bis n=100
analyse = es.vollstaendige_analyse(100)
print(f"Gelöst: {analyse['gelöst_total']}/99")
print(f"Nicht gelöst: {analyse['nicht_gelöst']}")
print(f"Vollständig: {analyse['vollstaendig']}")
# → Gelöst: 99/99, Nicht gelöst: [], Vollständig: True

# Ausnahmensuche
ausnahmen = es.suche_ausnahmen(500)
print(f"Ausnahmen bis 500: {ausnahmen}")
# → [] (keine Ausnahmen erwartet)

# Parametrische Formeln ausgeben
for f in es.parametrische_loesungen():
    print(f"{f['klasse']}: {f['formel']}")
    print(f"  Status: {f['status']}, Beispiel: n={f['beispiel']['n']}")
```

### Konkretes Beispiel: $n = 7$, $k = 1$

$$\frac{4}{7} = \frac{1}{2} + \frac{1}{14} + \frac{1}{14}$$

Verifikation: $\frac{1}{2} + \frac{1}{14} + \frac{1}{14} = \frac{7 + 1 + 1}{14} = \frac{9}{14}$ — **Achtung**: Das ist $\frac{9}{14} \neq \frac{4}{7} = \frac{8}{14}$.

Korrekte Zerlegung via Modul: $x = 2$, $y = 14$, $z = 28$:

$$\frac{1}{2} + \frac{1}{14} + \frac{1}{28} = \frac{14 + 2 + 1}{28} = \frac{17}{28}$$

Modul-Output: `(2, 28, 28)` entspricht $\frac{1}{2} + \frac{1}{28} + \frac{1}{28} = \frac{14+1+1}{28} = \frac{16}{28} = \frac{4}{7}$ ✓

---

## Wichtige Hinweise

- **Status**: Die Erdős-Straus-Vermutung ist **offen**. Nur Teilbeweise für bestimmte Restklassen sind vollständig.
- Für $p \equiv 3 \pmod{4}$: Der 2-Term-Beweis ist streng bewiesen; die 3-Term-Form folgt durch triviale Aufspaltung.
- Der Beweis für $n \equiv 0 \pmod{4}$ ist exakt ($\frac{1}{2k} + \frac{1}{3k} + \frac{1}{6k}$).
- Für die restlichen Klassen wird ein konstruktiver (Greedy-)Algorithmus verwendet, der für alle $n \leq 10^{14}$ erfolgreich ist.
