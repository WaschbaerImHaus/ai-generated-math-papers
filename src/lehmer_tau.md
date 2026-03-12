# lehmer_tau.py — Lehmer-Tau-Vermutung und Ramanujan-τ-Funktion

**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-12
**Modul:** `src/lehmer_tau.py`

---

## Inhaltsverzeichnis

1. [Mathematischer Hintergrund](#mathematischer-hintergrund)
2. [Klassen- und Methodenübersicht](#klassen--und-methodenübersicht)
3. [Wichtigste Algorithmen](#wichtigste-algorithmen)
4. [Beispielanwendungen](#beispielanwendungen)
5. [Kongruenztabelle](#kongruenztabelle)

---

## Mathematischer Hintergrund

### Die Ramanujan-τ-Funktion

Die **Ramanujan-Tau-Funktion** $\tau(n)$ ist definiert als der $n$-te Koeffizient der formalen Potenzreihenentwicklung der **Diskriminant-Modulform** $\Delta$:

$$\Delta(q) = q \cdot \prod_{n=1}^{\infty} (1 - q^n)^{24} = \sum_{n=1}^{\infty} \tau(n)\, q^n$$

$\Delta$ ist die eindeutige normalisierte **Spitzenform** (Cusp Form) vom Gewicht 12 bezüglich der vollen Modulgruppe $\mathrm{SL}_2(\mathbb{Z})$.

**Bekannte Werte:**

| $n$ | $\tau(n)$ |
|-----|-----------|
| 1 | 1 |
| 2 | -24 |
| 3 | 252 |
| 4 | -1472 |
| 5 | 4830 |
| 6 | -6048 |
| 7 | -16744 |
| 8 | 84480 |
| 9 | -113643 |
| 10 | -115920 |

### Lehmer-Tau-Vermutung

**Derrick Henry Lehmer (1947)** vermutete:

> $\tau(n) \neq 0$ für alle $n \geq 1$.

Diese Vermutung ist bis heute **offen**. Sie ist numerisch verifiziert für alle $n \leq 10^{22}$.

### Wichtige Eigenschaften

**Multiplikativität** (Hecke-Relationen):
$$\tau(mn) = \tau(m) \cdot \tau(n) \quad \text{falls } \gcd(m, n) = 1$$

Da $\tau$ multiplikativ ist, genügt es zum Beweis der Vermutung, $\tau(p) \neq 0$ für alle Primzahlen $p$ zu zeigen.

**Deligne-Schranke** (bewiesen 1974, Fields-Medaille 1978):
$$|\tau(p)| \leq 2 \cdot p^{11/2} \quad \text{für alle Primzahlen } p$$

Dies ist der berühmte Beweis der **Ramanujan-Vermutung** durch Pierre Deligne als Spezialfall der Weil-Vermutungen.

**Mod-691-Kongruenz** (fundamental):
$$\tau(n) \equiv \sigma_{11}(n) \pmod{691}$$

wobei $\sigma_{11}(n) = \sum_{d \mid n} d^{11}$. Die Zahl 691 taucht auf, weil sie den Zähler der Bernoulli-Zahl $B_{12} = -\tfrac{691}{2730}$ teilt.

### Notwendige Bedingungen für $\tau(p) = 0$

Nach Swinnerton-Dyer, Serre u.a. müssen folgende Kongruenzen **gleichzeitig** gelten, falls $\tau(p) = 0$:

$$\tau(p) \equiv 0 \pmod{691}$$
$$\tau(p) \equiv 0 \pmod{3}$$
$$\tau(p) \equiv 0 \pmod{5}$$
$$\tau(p) \equiv 0 \pmod{7}$$
$$\tau(p) \equiv 0 \pmod{23}$$

Bisher ist kein einziges $p$ bekannt, das alle diese Bedingungen erfüllt.

---

## Klassen- und Methodenübersicht

### Modul-Funktion (intern)

| Funktion | Signatur | Beschreibung |
|----------|----------|--------------|
| `_berechne_tau_liste` | `(n_max: int) → List[int]` | Berechnet $\tau(1), \ldots, \tau(n_{\max})$ via Polynommultiplikation |

### Klasse `LehmerTauAnalyse`

Hauptklasse zur Analyse der Lehmer-Tau-Vermutung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `berechne_tau` | `(n: int) → int` | Berechnet $\tau(n)$ exakt |
| `verifiziere_lehmer` | `(n_max: int) → Dict` | Sucht nach Nullstellen $\tau(n) = 0$ für $n \leq n_{\max}$ |
| `kongruenz_analyse` | `(p: int, mod_list?) → Dict` | Analysiert $\tau(p) \bmod m$ für verschiedene Moduli |
| `hecke_eigenvalue_schranke` | `(p: int) → Dict` | Prüft Deligne-Schranke $|\tau(p)| \leq 2p^{11/2}$ |
| `ausschluss_durch_kongruenzen` | `(p: int) → Dict` | Versucht $\tau(p) = 0$ auszuschließen |
| `bekannte_kongruenzen_table` | `() → Dict` | Tabelle aller bekannten Kongruenzbedingungen |

---

## Wichtigste Algorithmen

### 1. Koeffizientenextraktion via Polynommultiplikation

Die Berechnung von $\tau(n)$ nutzt die Produktdarstellung:

$$\Delta(q) = q \cdot \prod_{m=1}^{n_{\max}} (1 - q^m)^{24}$$

Dabei wird der Binomialsatz auf jedes Produkt angewendet:

$$(1 - q^m)^{24} = \sum_{k=0}^{24} \binom{24}{k} (-1)^k q^{km}$$

**Algorithmus:**
1. Initialisiere $p(q) = 1$ als Koeffizientenvektor der Länge $n_{\max}+1$
2. Für jedes $m = 1, \ldots, n_{\max}$: multipliziere $p(q)$ mit $(1-q^m)^{24}$
3. $\tau(n)$ ist der Koeffizient von $q^{n-1}$ in $p(q)$ (wegen des führenden $q$-Faktors)

Die Koeffizienten sind **exakte ganze Zahlen** — keine Rundungsfehler.

**Komplexität:** $O(n_{\max}^2)$ Multiplikationen

### 2. Verifikation der Lehmer-Vermutung

```python
ana = LehmerTauAnalyse()
erg = ana.verifiziere_lehmer(1000)
# {'verifiziert': True, 'nullstellen': [], 'schluss': '...'}
```

Da $\tau$ multiplikativ ist, kann man sich auf Primzahlpotenzen beschränken.

### 3. Kongruenzanalyse

Für jede Primzahl $p$ berechnet `kongruenz_analyse` die Reste $\tau(p) \bmod m$ für $m \in \{2, 3, 5, 7, 11, 13, 23, 691\}$ und prüft die bekannten notwendigen Bedingungen für $\tau(p) = 0$.

### 4. Ausschluss via mod-691-Kongruenz

Die Ramanujan-Kongruenz lautet:
$$\tau(p) \equiv 1 + p^{11} \pmod{691}$$

Falls $\tau(p) = 0$, muss gelten:
$$p^{11} \equiv -1 \equiv 690 \pmod{691}$$

---

## Beispielanwendungen

### τ(n) berechnen

```python
from lehmer_tau import LehmerTauAnalyse

ana = LehmerTauAnalyse()

# Einzelne Werte
print(ana.berechne_tau(1))   # 1
print(ana.berechne_tau(2))   # -24
print(ana.berechne_tau(5))   # 4830
```

### Lehmer-Vermutung verifizieren

```python
erg = ana.verifiziere_lehmer(500)
print(erg['verifiziert'])    # True
print(erg['nullstellen'])    # []
print(erg['schluss'])        # 'Lehmer-Vermutung verifiziert: τ(n) ≠ 0 für alle n ≤ 500.'
```

### Deligne-Schranke für Primzahlen

```python
for p in [2, 3, 5, 7, 11]:
    hs = ana.hecke_eigenvalue_schranke(p)
    print(f"p={p}: |τ(p)|={hs['abs_tau_p']}, Schranke≈{hs['deligne_schranke']:.0f}")
# p=2: |τ(p)|=24, Schranke≈91
# p=3: |τ(p)|=252, Schranke≈944
```

### Kongruenzanalyse für p=691

```python
k = ana.kongruenz_analyse(691)
print(k['tau_p'])                    # τ(691)
print(k['könnte_null_sein'])         # False (für bekannte Primzahlen)
print(k['kongruenzen'][691])         # {'tau_mod_m': ..., 'null': False}
```

### Mod-691-Tabelle verifizieren

```python
tab = ana.bekannte_kongruenzen_table()
for v in tab['verifikation_erste_n'][:5]:
    print(f"n={v['n']}: τ={v['tau_n']}, σ₁₁={v['sigma_11_n']}, "
          f"mod-691 stimmt: {v['mod_691_stimmt']}")
```

---

## Kongruenztabelle

| Modul | Kongruenz | Quelle | Notwendig für $\tau(p)=0$ |
|-------|-----------|--------|--------------------------|
| 2 | $\tau(n) \equiv \sigma_{11}(n) \pmod{8}$ | Ramanujan 1916 | $\tau \equiv 0 \pmod 8$ |
| 3 | $\tau(n) \equiv \sigma_{11}(n) \pmod{3}$ | Ramanujan 1916 | $\sigma_{11}(n) \equiv 0 \pmod 3$ |
| 5 | $\tau(n) \equiv \sigma_{11}(n) \pmod{5}$ | Ramanujan 1916 | $\sigma_{11}(n) \equiv 0 \pmod 5$ |
| 7 | $\tau(n) \equiv \sigma_{11}(n) \pmod{7}$ | Ramanujan 1916 | $\sigma_{11}(n) \equiv 0 \pmod 7$ |
| 23 | $\tau(p) \equiv 0$ oder $\pm 2 \pmod{23}$ | Swinnerton-Dyer 1973 | $\tau(p) \equiv 0 \pmod{23}$ |
| 691 | $\tau(n) \equiv \sigma_{11}(n) \pmod{691}$ | Ramanujan 1916 | $\sigma_{11}(n) \equiv 0 \pmod{691}$ |

**Erklärung der 691:** Die Bernoulli-Zahl $B_{12} = -\frac{691}{2730}$ hat 691 als Zähler. Die Eisenstein-Reihe $E_{12}$ und $\Delta$ stimmen modulo 691 überein, was die Kongruenz erklärt.

---

## Verbindung zur Modulform-Theorie

$\Delta$ liegt im Raum $S_{12}(\mathrm{SL}_2(\mathbb{Z}))$ der Spitzenformen vom Gewicht 12. Da $\dim S_{12} = 1$, ist $\Delta$ die **eindeutige** normalisierte Spitzenform dieses Gewichts. Die Hecke-Operatoren $T_p$ wirken auf $\Delta$ durch Skalierung mit $\tau(p)$:

$$T_p \Delta = \tau(p) \cdot \Delta$$

Die Lehmer-Vermutung besagt also, dass $\Delta$ kein **Hecke-Eigenvektor zum Eigenwert 0** ist — was aus rein algebraischen Gründen möglich wäre, aber offenbar nie eintritt.

---

*Dokumentation generiert für das specialist-maths Projekt. Autor: Michael Fuhrmann.*
