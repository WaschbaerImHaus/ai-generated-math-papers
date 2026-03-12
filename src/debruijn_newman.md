# debruijn_newman.py — de Bruijn-Newman-Konstante Λ und Riemann-Hypothese

**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-12
**Modul:** `src/debruijn_newman.py`

---

## Inhaltsverzeichnis

1. [Mathematischer Hintergrund](#mathematischer-hintergrund)
2. [Klassen- und Methodenübersicht](#klassen--und-methodenübersicht)
3. [Wichtigste Algorithmen](#wichtigste-algorithmen)
4. [Beispielanwendungen](#beispielanwendungen)
5. [Historische Schranken](#historische-schranken)

---

## Mathematischer Hintergrund

### Die Kernfunktion Φ(u)

Die Theorie basiert auf der Funktion:

$$\Phi(u) = \sum_{n=1}^{\infty} \left(2\pi^2 n^4 e^{9u} - 3\pi n^2 e^{5u}\right) \exp\!\left(-\pi n^2 e^{4u}\right)$$

Diese Funktion entsteht aus der **Jacobi-Theta-Funktion**
$$\omega(x) = \sum_{n=1}^{\infty} e^{-\pi n^2 x}$$
durch zweifache Differentiation und eine geeignete Substitution. Sie fällt für $u \to \infty$ superexponentiell ab (dominiert durch $e^{-\pi e^{4u}}$), sodass nur wenige Terme ($n=1,2$) praktisch relevant sind.

### Die parametrische Funktionsfamilie $H_t$

de Bruijn (1950) betrachtete:

$$H_t(x) = \int_0^{\infty} e^{t u^2} \cdot \Phi(u) \cdot \cos(xu)\, du \quad (t \in \mathbb{R},\ x \in \mathbb{R})$$

Für $t = 0$ ergibt sich ein besonderer Spezialfall:

$$H_0(x) = \frac{1}{2}\, \xi\!\left(\tfrac{1}{2} + ix\right)$$

wobei $\xi(s) = \frac{1}{2} s(s-1) \pi^{-s/2} \Gamma(s/2) \zeta(s)$ die **vollständige Riemann-Zeta-Funktion** ist.

Die Nullstellen von $H_0$ entsprechen genau den **imaginären Teilen** der nicht-trivialen Nullstellen der Riemann-Zeta-Funktion auf der kritischen Geraden $\mathrm{Re}(s) = \frac{1}{2}$:

$$y_1 \approx 14.134725,\quad y_2 \approx 21.022040,\quad y_3 \approx 25.010857,\quad \ldots$$

### Die de Bruijn-Newman-Konstante Λ

**de Bruijn (1950)** bewies:
- Wenn $H_{t_0}$ nur reelle Nullstellen hat, dann hat $H_t$ nur reelle Nullstellen für alle $t \geq t_0$.
- Damit ist $\Lambda$ wohldefiniert:

$$\Lambda = \inf\{ t \in \mathbb{R} : H_t \text{ hat nur reelle Nullstellen}\}$$

### Verbindung zur Riemann-Hypothese

Die zentrale Äquivalenz lautet:

$$\textbf{RH} \iff \Lambda \leq 0$$

**Vollständige Äquivalenzkette:**
1. RH: Alle nicht-trivialen Nullstellen von $\zeta(s)$ liegen auf $\mathrm{Re}(s) = \frac{1}{2}$
2. $\iff$ $\xi(\frac{1}{2} + ix)$ hat nur reelle Nullstellen
3. $\iff$ $H_0(x) = \frac{1}{2}\xi(\frac{1}{2}+ix)$ hat nur reelle Nullstellen
4. $\iff$ $\Lambda \leq 0$

### Newman-Vermutung und ihr Beweis

**C.M. Newman (1976)** vermutete:
$$\Lambda \geq 0$$

**Bedeutung:** Falls RH wahr ist, dann gilt $\Lambda = 0$ (nicht $\Lambda < 0$). Die RH wäre eine „scharfe" Aussage.

**Rodgers & Tao (2018)** bewiesen die Newman-Vermutung (arXiv:1801.05914):
$$\Lambda \geq 0 \quad \textbf{BEWIESEN}$$

**Aktueller Wissensstand (2026):**
$$0 \leq \Lambda \leq 0.2$$

Die untere Schranke stammt von Rodgers-Tao (2018), die obere Schranke von Platt & Trudgian (2021).

---

## Klassen- und Methodenübersicht

### Modul-Konstante

```python
BEKANNTE_NULLSTELLEN_GAMMA: List[float]  # Erste 10 Riemann-Nullstellen (Im-Teile)
```

### Klasse `DeBruijnNewman`

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `__init__` | `(praezision=50)` | Setzt mpmath-Präzision |
| `phi_funktion` | `(u, N_max=20) → float` | Berechnet $\Phi(u)$ |
| `H_t` | `(x, t, N_max=30) → float` | Berechnet $H_t(x)$ |
| `suche_nullstellen_H_t` | `(t, x_bereich?, N_max?, schritte?) → List[float]` | Nullstellen von $H_t$ via Bisektion |
| `ist_alle_reell` | `(t, grenze?, N_max?, toleranz?) → bool` | Heuristische Reellwertprüfung |
| `schranke_lambda_oben` | `(kandidat_t) → Dict` | Numerische Evidenz für $\Lambda \leq$ kandidat_t |
| `riemann_zusammenhang` | `() → Dict` | Formale Beschreibung der RH-Äquivalenz |
| `newman_vermutung_status` | `() → Dict` | Status-Bericht zur Newman-Vermutung |
| `bekannte_schranken_history` | `() → List[Dict]` | Historische Schranken-Entwicklung |
| `nullstellen_H0_berechnen` | `(x_max?, N_max?) → List[float]` | Erste Nullstellen von $H_0$ numerisch |
| `vergleiche_mit_riemann_nullstellen` | `(nullstellen) → Dict` | Vergleich mit bekannten Werten |
| `visualisiere_H0` | `(x_max?, N_max?, dateiname?) → bool` | matplotlib-Plot von $H_0$ |

**Interne Hilfsmethoden:**

| Methode | Beschreibung |
|---------|--------------|
| `_phi_mpmath` | Hochpräzise $\Phi(u)$-Berechnung via mpmath |
| `_phi_python` | Python-Fallback für $\Phi(u)$ |
| `_H_t_mpmath` | Hochpräzise $H_t$-Berechnung via `mpmath.quad` |
| `_H_t_quadratur` | Simpson-Quadratur als Python-Fallback |
| `_H_t_quadratur_mp` | mpmath-Simpson als zweiter Fallback |

---

## Wichtigste Algorithmen

### 1. Berechnung von Φ(u)

Für jedes $n = 1, 2, \ldots, N_{\max}$:
1. Berechne den Dämpfungsfaktor $D_n = \exp(-\pi n^2 e^{4u})$
2. Falls $D_n < 10^{-\text{präzision}+5}$: abbrechen (Rest vernachlässigbar)
3. Berechne Wachstumsterm $W_n = 2\pi^2 n^4 e^{9u} - 3\pi n^2 e^{5u}$
4. Addiere $W_n \cdot D_n$ zum Ergebnis

Für $u \geq 0$ dominiert der Term $n=1$ vollständig; $n=2$ trägt nur für kleine $u$ bei.

### 2. Numerische Integration von H_t

**Methode 1 (bevorzugt): mpmath adaptive Quadratur**

```python
H_t(x) = mpmath.quad(lambda u: exp(t*u**2) * Phi(u) * cos(x*u), [0, 5])
```

Das Integrationsintervall $[0, 5]$ reicht aus, da $\Phi(u) \approx 0$ für $u > 2$.

**Methode 2 (Fallback): Simpson-Quadratur**

Mit $M = 500$ Intervallen der Breite $h = 5/500$:

$$H_t(x) \approx \frac{h}{3} \sum_{i=0}^{M} w_i \cdot e^{t(ih)^2} \cdot \Phi(ih) \cdot \cos(x \cdot ih)$$

wobei $w_i \in \{1, 4, 2\}$ die Simpson-Gewichte sind.

### 3. Nullstellensuche via Vorzeichenwechsel + Bisektion

```
Bewerte H_t an schritte Gitterpunkten in [a,b]
Für jedes Intervall mit f_prev * f_curr < 0:
    Bisektion (50 Schritte) → Präzision ~15 Dezimalstellen
```

Für $t = 0$ sollten die Nullstellen die bekannten Riemann-Nullstellen reproduzieren.

### 4. Heuristische Reell-Prüfung

Für $t \approx 0$: Vergleiche gefundene Nullstellen mit den bekannten Werten $\gamma_1, \gamma_2, \ldots$ Abweichung $\leq 1.0$ wird akzeptiert (numerische Näherung).

---

## Beispielanwendungen

### RH-Äquivalenz abrufen

```python
from debruijn_newman import DeBruijnNewman

dbn = DeBruijnNewman(praezision=30)

zusammenhang = dbn.riemann_zusammenhang()
print(zusammenhang['hauptsatz'])
# "RH (Riemann-Hypothese) ⟺ Λ ≤ 0"
print(zusammenhang['rh_aequivalenz'])
```

### Newman-Vermutung-Status

```python
status = dbn.newman_vermutung_status()
print(status['status'])              # "BEWIESEN ✓"
print(status['beweis_quelle'])       # "B. Rodgers & T. Tao (2018)"
print(status['aktueller_wissenstand']['obere_schranke'])
# "Λ ≤ 0.2  (Platt-Trudgian 2021, rigoros)"
```

### Phi-Funktion auswerten

```python
for u in [0.0, 0.5, 1.0, 2.0]:
    phi = dbn.phi_funktion(u, N_max=10)
    print(f"Phi({u}) = {phi:.6e}")
# Phi(0.0) ≈ 2.5e-1  (dominiert von n=1-Term)
# Phi(2.0) ≈ 0       (superexponentiell abfallend)
```

### H_0(x) auswerten

```python
# H_0 ist proportional zu xi(1/2 + ix)
# Hat Nullstellen bei den Riemann-Nullstellen
h0_bei_0 = dbn.H_t(0.0, t=0.0)
h0_bei_y1 = dbn.H_t(14.134725, t=0.0)
print(f"H_0(0)    = {h0_bei_0:.6f}")
print(f"H_0(y_1)  = {h0_bei_y1:.6f}")  # nahe 0
```

### Nullstellen von H_0 berechnen

```python
nullstellen = dbn.nullstellen_H0_berechnen(x_max=50.0, N_max=15)
print(f"Gefundene Nullstellen: {len(nullstellen)}")
bekannt = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]
for i, y in enumerate(nullstellen[:5]):
    print(f"  y_{i+1}: gefunden={y:.6f}  bekannt={bekannt[i]:.6f}")
```

### Mit bekannten Riemann-Nullstellen vergleichen

```python
vergleich = dbn.vergleiche_mit_riemann_nullstellen(nullstellen)
print(vergleich['zusammenfassung'])
for r in vergleich['ergebnisse'][:3]:
    print(f"  γ_{r['index']}: erwartet={r['erwartet']:.6f}, "
          f"gefunden={r['gefunden']:.6f}, Fehler={r['absoluter_fehler']:.6f}")
```

### Obere Schranke numerisch prüfen

```python
for t_test in [0.0, 0.1, 0.2, 0.5]:
    schranke = dbn.schranke_lambda_oben(t_test)
    print(f"t={t_test}: {schranke['schlussfolgerung']}")
```

### Historische Schranken ausgeben

```python
history = dbn.bekannte_schranken_history()
for e in history:
    print(f"{e['jahr']}: {e['schranke']}")
```

### H_0 visualisieren

```python
# Zeichnet H_0(x) mit markierten Riemann-Nullstellen
dbn.visualisiere_H0(x_max=50.0, dateiname="H0_plot.png")
```

---

## Historische Schranken

| Jahr | Schranke | Richtung | Quelle |
|------|----------|----------|--------|
| 1950 | $\Lambda \leq \frac{1}{2}$ | oben | N.G. de Bruijn |
| 1976 | $\Lambda \geq -\infty$ (Vermutung: $\Lambda \geq 0$) | unten | C.M. Newman |
| 1988 | $\Lambda < \frac{1}{2}$ (scharf) | oben | Newman & Goldberg |
| 2009 | $\Lambda \geq -2.7 \cdot 10^{-9}$ | unten | Saouter, Gourdon, Demichel |
| **2018** | $\Lambda \geq 0$ **BEWIESEN** | unten | **Rodgers & Tao (arXiv:1801.05914)** |
| 2020 | $0 \leq \Lambda \leq 0.22$ | kombiniert | Polymath15 (T. Tao et al.) |
| 2021 | $0 \leq \Lambda \leq 0.2$ | kombiniert | Platt & Trudgian |

---

## Mathematische Tiefe: Warum gilt RH ⟺ Λ ≤ 0?

**de Bruijn's Schlüsselresultat:** Der Operator "Faltung mit Gauß-Kern $e^{t u^2}$" verschiebt die Nullstellen von $H$ in die komplexe Ebene, sobald $t < \Lambda$. Genauer: $H_{t_0}$ hat nur reelle Nullstellen $\implies H_t$ hat nur reelle Nullstellen für alle $t \geq t_0$. Der "Schwellenwert" $\Lambda$ ist das Infimum dieser $t$-Werte.

**Spezialfall $t = 0$:** $H_0 = \frac{1}{2}\xi(\frac{1}{2} + ix)$. Die RH besagt genau, dass $\xi(\frac{1}{2}+ix)$ nur reelle Nullstellen hat, d.h. $\Lambda \leq 0$.

**Rodgers-Tao 2018:** Zeigten, dass für $t < 0$ immer nicht-reelle Nullstellen entstehen (probabilistisch-analytischer Beweis über GUE-Korrelationen der Nullstellenverteilung). Damit ist $\Lambda \geq 0$ bewiesen.

**Konsequenz:** Die RH ist äquivalent zu $\Lambda = 0$ (denn $\Lambda \geq 0$ ist bewiesen und $\Lambda \leq 0$ wäre äquivalent zur RH). Die aktuelle Schranke $\Lambda \leq 0.2$ zeigt, dass wir der RH aus dieser Richtung sehr nahe kommen.

---

## Abhängigkeiten

| Bibliothek | Verwendung | Optional |
|------------|------------|----------|
| `mpmath` | Hochpräzise Arithmetik, adaptive Quadratur | Nein (empfohlen) |
| `numpy` | Vektoroperationen | Ja (aktuell ungenutzt) |
| `matplotlib` | Visualisierung von $H_0$ | Ja |

---

*Dokumentation generiert für das specialist-maths Projekt. Autor: Michael Fuhrmann.*
