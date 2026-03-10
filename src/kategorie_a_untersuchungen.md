# kategorie_a_untersuchungen.py – Dokumentation

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Build:** 38

---

## Überblick

Dieses Modul implementiert formale mathematische Untersuchungen zu den
**Kategorie-A-Vermutungen** – also jenen offenen Problemen, für die ein
Beweis in absehbarer Zeit als wahrscheinlich gilt, weil:

- partielle Beweise oder enge Näherungen existieren,
- bedingte Beweise (unter GRH o.ä.) vorliegen,
- numerische Verifikationen bis sehr großen Grenzen durchgeführt wurden.

---

## Struktur

### Hilfsfunktionen

| Funktion | Beschreibung |
|----------|-------------|
| `_sieve_primes(limit)` | Sieb des Eratosthenes; gibt alle Primzahlen ≤ `limit` zurück |
| `_is_prime_fast(n)` | Schneller Primzahltest via SymPy (Miller-Rabin, deterministisch) |

---

## Klassen

### 1. `GoldbachUntersuchung`

**Vermutung (1742):** Jede gerade Zahl $n > 2$ ist Summe zweier Primzahlen:
$$n = p + q, \quad p, q \in \mathbb{P}$$

**Bekannte Resultate:**
- Helfgott (2013): Ternäre Goldbach vollständig bewiesen
- Chen (1973): $n = p + P_2$ für alle hinreichend großen geraden $n$
- Numerisch verifiziert bis $4 \times 10^{18}$

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `verifiziere_bis(grenze)` | Verifiziert die Goldbach-Vermutung für alle geraden Zahlen bis `grenze` |
| `goldbach_zerlegungen_anzahl(n)` | Zählt alle Primzahlpaare $p + q = n$ mit $p \leq q$ |
| `hardy_littlewood_singular_series(n)` | Berechnet $G(n) \approx 2C_2 \prod_{p\|n,p>2} \frac{p-1}{p-2} \cdot \frac{n}{(\ln n)^2}$ |
| `vinogradov_schranke()` | Beschreibung von Vinogradovs bewiesenem Dreiprimzahlsatz |
| `chen_satz()` | Beschreibung von Chens bewiesenem Satz (1973) |

---

### 2. `ZwillingsprimzahlUntersuchung`

**Vermutung (~1849):** Es gibt unendlich viele Primzahlpaare $(p, p+2)$.

**Bekannte Resultate:**
- Brun (1919): $B_2 = \sum_{(p,p+2)\text{ Zwilling}} (\tfrac{1}{p} + \tfrac{1}{p+2}) < \infty$
- Zhang (2013): $\liminf_{n \to \infty}(p_{n+1}-p_n) < 7 \times 10^7$
- Maynard (2013): $< 600$; Polymath8b (2014): $< 246$

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `verifiziere_bis(grenze)` | Findet alle Zwillingsprimzahlpaare $(p, p+2)$ bis `grenze` |
| `brun_konstante_naerung(grenze)` | Numerische Näherung der Brun-Konstante $B_2 \approx 1.902$ |
| `zhang_maynard_ergebnis()` | Zusammenfassung der Durchbrüche von Zhang und Maynard |
| `dichte_zwillingsprimzahlen(x)` | Hardy-Littlewood-Schätzung: $\pi_2(x) \sim 2C_2 \int_2^x \frac{dt}{(\ln t)^2}$ |

---

### 3. `LegendreUntersuchung`

**Vermutung (1798):** Für jedes $n \geq 1$ liegt eine Primzahl im Intervall $(n^2, (n+1)^2)$.

**Bekannte Resultate:**
- Bertrand-Postulat (Tschebyschow 1850): **BEWIESEN** – Zwischen $n$ und $2n$ liegt stets eine Primzahl
- Ingham (1937): Primzahl zwischen $n^3$ und $(n+1)^3$ (bewiesen)
- Huxley (1972): $\pi(x + x^\theta) - \pi(x) \sim x^\theta / \ln x$ für $\theta > 7/12$

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `verifiziere_bis(n_max)` | Verifiziert Legendres Vermutung für $n = 1, \ldots, n_{\max}$ |
| `bertrand_postulat(n)` | Findet Primzahl $p$ mit $n < p \leq 2n$ (Tschebyschow, bewiesen) |
| `huxley_schranke(x, theta)` | Schätzt $\pi(x+x^\theta) - \pi(x)$ nach Huxley |
| `ingham_kubisch(n)` | Findet Primzahl zwischen $n^3$ und $(n+1)^3$ (Ingham, bewiesen) |

---

### 4. `BrocardUntersuchung`

**Vermutung (1904):** Für alle $n \geq 2$ gilt:
$$\pi(p_{n+1}^2) - \pi(p_n^2) \geq 4$$
d.h. zwischen den Quadraten zweier aufeinanderfolgender Primzahlen liegen mindestens 4 Primzahlen.

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `verifiziere_bis(primzahl_index_max)` | Verifiziert Brocards Vermutung bis zum $n$-ten Primzahl-Index |
| `verbindung_zu_legendre()` | Erklärt die Verbindung zur Legendres Vermutung |

---

### 5. `ErdosStrausUntersuchung`

**Vermutung (Erdős, 1948):** Für jedes $n \geq 2$ existieren positive ganze Zahlen $x, y, z$ mit:
$$\frac{4}{n} = \frac{1}{x} + \frac{1}{y} + \frac{1}{z}$$

**Bekannte Resultate:**
- Verifiziert für alle $n \leq 10^{14}$ (Elsholtz-Tao, 2013)
- Partieller Beweis via Restklassen (für $n \equiv 0 \pmod{4}$ vollständig)

**Schlüssel-Identität** für $n \equiv 0 \pmod{4}$, d.h. $n = 4k$:
$$\frac{4}{4k} = \frac{1}{k} = \frac{1}{2k} + \frac{1}{3k} + \frac{1}{6k}$$

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `finde_zerlegung(n)` | Sucht Zerlegung $4/n = 1/x + 1/y + 1/z$ (gierig + brute force) |
| `_brute_force(n)` | Erschöpfende Suche für $n \leq 10000$ |
| `verifiziere_bis(grenze)` | Verifiziert die Vermutung für $n = 2, \ldots, \text{grenze}$ |
| `beweis_via_restklassen()` | Partieller Beweis für Restklassen mod 4 |

---

### 6. `ArtinVermutungUntersuchung`

**Vermutung (Artin, 1927):** Sei $a \in \mathbb{Z}$ mit $a \neq 0, \pm 1$ und $a$ kein vollständiges Quadrat. Dann ist $a$ primitive Wurzel modulo $p$ für unendlich viele Primzahlen $p$, mit Dichte:
$$A(a) = \prod_p \left(1 - \frac{1}{p(p-1)}\right) \approx 0.3739558...$$

**Beweis unter GRH** (Hooley, 1967): Vollständig bedingt bewiesen.
**Ohne GRH** (Gupta-Murty-Heath-Brown, 1984): Gilt für alle $a$ außer höchstens zwei Ausnahmen.

**Klassenattribut:** `ARTIN_KONSTANTE = 0.3739558136192022`

#### Methoden

| Methode | Beschreibung |
|---------|-------------|
| `ist_primitive_wurzel(a, p)` | Prüft ob $\text{ord}_p(a) = p-1$ (d.h. $a$ ist primitive Wurzel mod $p$) |
| `dichte_verifizieren(a, bis_primzahl)` | Beobachtete vs. vorhergesagte Dichte der primitiven Wurzeln |
| `artin_konstante_berechnen(n_primes)` | Numerische Näherung von $A = \prod_p (1 - 1/(p(p-1)))$ |
| `hooley_bedingung()` | Beschreibung von Hooleys bedingtem Beweis (unter GRH) |

---

### 7. `GoldbachPartiellerBeweis`

Beweist Goldbach für spezielle Zahlenklassen vollständig.

| Methode | Beschreibung |
|---------|-------------|
| `goldbach_fuer_vielfache_von_6(k)` | Goldbach für $n = 6k$: Findet Paar $(p, q)$ mit $p + q = 6k$ |
| `goldbach_fuer_zweierpotenzen(k)` | Goldbach für $n = 2^k$: Findet Paar $(p, q)$ mit $p + q = 2^k$ |
| `satz_gerade_zahlen_als_summe_von_zwei_halbprimzahlen(n)` | Findet alle $P_2$-Paare $(p, q)$ mit $p + q = n$ |

---

## Modulare Zusammenfassung

```python
from kategorie_a_untersuchungen import kategorie_a_zusammenfassung

summary = kategorie_a_zusammenfassung()
# summary['Goldbach']['status'] → 'Offen'
# summary['Artin']['status'] → 'BEWIESEN unter GRH ...'
```

---

## Mathematische Zusammenhänge

```
Goldbach (stark)  →  Goldbach (schwach, ternär)
                              ↑ Helfgott 2013 BEWIESEN

Legendre          →  Brocard (quantitativ stärker)
Legendre     ←  (bedingt, unter RH)
Bertrand          ←  Legendre (Spezialfall)
Bertrand          ↑ Tschebyschow 1850 BEWIESEN

Artin              ↑ Hooley 1967 BEWIESEN unter GRH
                    ↑ GMH 1984: für fast alle a UNBEDINGT
```

---

## Beispielnutzung

```python
from kategorie_a_untersuchungen import (
    GoldbachUntersuchung, ZwillingsprimzahlUntersuchung,
    ErdosStrausUntersuchung, ArtinVermutungUntersuchung
)

# Goldbach bis 1000 verifizieren
gb = GoldbachUntersuchung()
result = gb.verifiziere_bis(1000)
print(result['verifiziert'])  # True

# Zwillingsprimzahlpaare bis 100
zp = ZwillingsprimzahlUntersuchung()
print(zp.verifiziere_bis(100)['erste_10'])

# 4/97 = 1/x + 1/y + 1/z
es = ErdosStrausUntersuchung()
x, y, z = es.finde_zerlegung(97)
print(f"4/97 = 1/{x} + 1/{y} + 1/{z}")

# Artin-Konstante berechnen
artin = ArtinVermutungUntersuchung()
print(artin.artin_konstante_berechnen(200))  # ≈ 0.3740
```

---

## Tests

Die Testdatei `tests/test_kategorie_a_untersuchungen.py` enthält **102 Tests**
in 10 Testklassen, die alle relevanten Methoden und Edge-Cases abdecken.
