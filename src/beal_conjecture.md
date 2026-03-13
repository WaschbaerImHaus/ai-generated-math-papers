# beal_conjecture.py — Erläuterungen

**Autor:** Michael Fuhrmann
**Datum:** 2026-03-13
**Status:** CONJECTURE (offen seit 1993, Preisgeld: 1 Mio. USD)

---

## Überblick

Dieses Modul implementiert die computationale Verifikation von **Beal's Vermutung** und analysiert ihre Verbindung zur **ABC-Vermutung**.

### Die Vermutung (CONJECTURE)

> Falls $A^x + B^y = C^z$ mit $A, B, C, x, y, z \in \mathbb{Z}^+$ und $x, y, z \geq 3$, dann haben $A$, $B$, $C$ einen gemeinsamen Primteiler.

Äquivalent: $\gcd(A, B, C) > 1$.

**Preisgeld:** 1.000.000 USD (Andrew Beal, Stand 2023)
**Status:** Offen. Weder Beweis noch Gegenbeispiel bekannt.

---

## Mathematischer Hintergrund

### Verbindung zu Fermat's Letztem Satz

Fermat's Letzter Satz (Wiles 1995) ist ein Spezialfall: $x = y = z = n \geq 3$, $\gcd(A,B,C) = 1$ kann nicht gelten, da keine Lösung existiert. Beal verallgemeinert auf unterschiedliche Exponenten.

### Bekannte Familien gültiger Tripel

Alle gültigen Beal-Tripel haben gemeinsamen Primteiler:

| Tripel | Wert | Gemeinsamer Faktor |
|--------|------|--------------------|
| $2^3 + 2^3 = 2^4$ | $8+8=16$ | 2 |
| $2^6 + 2^6 = 2^7$ | $64+64=128$ | 2 |
| $3^3 + 6^3 = 3^5$ | $27+216=243$ | 3 |
| $27^3 + 84^3 = 3^{11}$ | $\ldots$ | 3 |

### Verbindung zur ABC-Vermutung

Die **ABC-Vermutung** (Oesterlé-Masser 1985, CONJECTURE) besagt:
Für jedes $\varepsilon > 0$ existiert $K_\varepsilon > 0$, sodass für alle teilerfremden $a + b = c$:
$$c < K_\varepsilon \cdot \text{rad}(abc)^{1+\varepsilon}$$

Die **Qualität** eines ABC-Tripels ist:
$$q(a,b,c) = \frac{\log c}{\log \text{rad}(abc)}$$

Falls Beal falsch wäre (teilerfremdes $A^x + B^y = C^z$), wäre:
$$q = \frac{z \log C}{\log \text{rad}(A \cdot B \cdot C)}$$

Für $x, y, z \geq 3$ und $\gcd(A,B,C) = 1$ erhält man hohe Qualitätswerte, was ABC einschränkt.

---

## Implementierte Klassen

### `BealTriple`

Repräsentiert $A^x + B^y = C^z$.

| Methode | Beschreibung |
|---------|--------------|
| `is_valid()` | Prüft $A^x + B^y = C^z$ arithmetisch |
| `has_common_factor()` | $\gcd(A,B,C) > 1$? |
| `satisfies_beal()` | Beal für dieses Tripel erfüllt? |
| `radical()` | $\text{rad}(A \cdot B \cdot C)$ |
| `common_factor()` | $\gcd(A,B,C)$ |

### `BealChecker`

Effiziente Suche via Lookup-Tabellen.

**Algorithmus:**
1. Berechne alle $k^e$ für $k \leq N$, $e \in \{3, \ldots, \text{max\_exp}\}$
2. Speichere in Dict: Wert $\to$ [(Basis, Exp), ...]
3. Für jedes Paar $(A^x, B^y)$: prüfe ob $A^x + B^y$ im Dict als $C^z$ vorkommt
4. Sammle alle gültigen Tripel

| Methode | Beschreibung |
|---------|--------------|
| `precompute_powers()` | Baut Lookup-Tabelle auf |
| `find_all_triples()` | Alle Beal-Tripel mit Basen $\leq N$ |
| `find_counterexamples()` | Gegenbeispiele (erwartet: leer) |
| `known_families()` | Vordefinierte gültige Familien |

### `ABCConnection`

Analyse der ABC-Qualitäten.

| Methode | Beschreibung |
|---------|--------------|
| `radical(n)` | $\text{rad}(n) = \prod_{p \mid n} p$ |
| `quality(a,b,c)` | $q = \log c / \log \text{rad}(abc)$ |
| `beal_abc_bound(A,x,B,y,C,z)` | ABC-Qualität für Beal-Tripel |
| `high_quality_abc_triples(limit)` | Tripel mit $q > 1$ |

---

## Das Radikal

$$\text{rad}(n) = \prod_{\substack{p \mid n \\ p \text{ prim}}} p$$

Beispiele:
- $\text{rad}(12) = \text{rad}(2^2 \cdot 3) = 2 \cdot 3 = 6$
- $\text{rad}(p^k) = p$ für Primzahl $p$
- $\text{rad}(A^x \cdot B^y \cdot C^z) = \text{rad}(A \cdot B \cdot C)$

---

## Suchschranken

Für $N = 100$, $\text{max\_exp} = 10$:
- Potenzen: $100 \times 8 = 800$ Einträge
- Paare: $O(800^2)$ Kombinationen
- Gefundene Tripel: mehrere hundert

**Ergebnis bisher (numerisch):** Kein Gegenbeispiel für $A, B, C \leq 1000$ bekannt.

---

## Literatur

- Beal, A.R. (1993): Erste Formulierung der Vermutung
- Granville, A. (2002): Verbindung zu ABC, *Bull. AMS*
- Mauldin, R.D. (1997): *A generalization of Fermat's Last Theorem*, *Notices AMS*
- Wiles, A. (1995): Beweis von FLT (Spezialfall von Beal), *Ann. Math.*
- Mochizuki, S. (2012): Inter-universal Teichmüller theory (Beweis-Versuch für ABC, disputed)
