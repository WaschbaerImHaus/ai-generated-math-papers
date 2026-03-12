# riemann_siegel_ext.py — Erweiterte Riemann-Siegel-Formel

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Datum**: 2026-03-12

---

## Überblick

Dieses Modul implementiert die **Riemann-Siegel-Formel** zur effizienten Berechnung
der Riemann-Zeta-Funktion auf der kritischen Geraden Re(s) = 1/2.

Die zentrale Größe ist die **Hardy-Z-Funktion**:

$$Z(t) = e^{i\theta(t)} \cdot \zeta(1/2 + it)$$

Da $e^{i\theta(t)}$ auf der kritischen Geraden die Phase von $\zeta$ kompensiert,
ist $Z(t)$ für reelles $t$ stets **reell**. Nullstellen von $Z(t)$ entsprechen
Nullstellen von $\zeta(s)$ auf der kritischen Geraden.

---

## Mathematischer Hintergrund

### Siegelsche θ-Funktion

$$\theta(t) = \operatorname{Im}\bigl(\ln \Gamma(1/4 + it/2)\bigr) - \frac{t}{2} \ln \pi$$

- $\theta(t)$ ist monoton wachsend für $t > t_{\min} \approx 6.29$
- $\theta(t) > 0$ erst für $t > g_0 \approx 17.845$
- **Stirling-Asymptotik** für große $t$:
  $$\theta(t) \approx \frac{t}{2}\ln\frac{t}{2\pi} - \frac{t}{2} - \frac{\pi}{8} + \frac{1}{48t} + O(t^{-3})$$

### Riemann-Siegel-Hauptformel

Sei $N = \lfloor\sqrt{t/(2\pi)}\rfloor$. Dann gilt:

$$Z(t) = 2 \sum_{n=1}^{N} \frac{\cos(\theta(t) - t \ln n)}{\sqrt{n}} + R(t)$$

**Fehlerterm** (Gabcke 1979):
$$R(t) = O(t^{-1/4}), \quad |R(t)| \leq 0.053 \cdot t^{-3/4} + 0.14 \cdot t^{-5/4}$$

Der Hauptterm $C_0$ des Fehlerterms ist:
$$R(t) \approx (-1)^{N-1} \cdot t^{-1/4} \cdot C_0(p)$$
mit $p = \{\sqrt{t/(2\pi)}\}$ (gebrochener Anteil) und:
$$C_0(p) = \frac{\cos(2\pi(p^2 - p - 1/16))}{\cos(2\pi p)}$$

---

## Gram-Punkte und Gram-Gesetz

### Definition

Der $n$-te **Gram-Punkt** $g_n$ ist definiert durch:
$$\theta(g_n) = n\pi$$

Bekannte Werte: $g_0 \approx 17.845$, $g_1 \approx 23.170$, $g_2 \approx 27.670$, ...

### Gram-Gesetz

**Gram's Law** (J.-P. Gram, 1903): $Z(g_n)$ sollte das Vorzeichen $(-1)^n$ haben.

- $n$ gerade: $Z(g_n) > 0$ erwartet
- $n$ ungerade: $Z(g_n) < 0$ erwartet

**Gilt in ~73% aller Fälle** (Titchmarsh 1935). Die verbleibenden ~27% sind
**Gram-Ausnahmen** ("Gram failures").

### Gram-Blöcke und Rosser-Regel

Wenn Gram-Intervalle die Gram-Law-Bedingung nicht einzeln erfüllen, werden sie
zu **Gram-Blöcken** zusammengefasst. Die **Rosser-Regel** (1941) besagt:

> Wenn $Z(g_n) > 0$ (für gerades $n$) und $Z(g_{n+2}) > 0$, dann liegt
> **genau eine** Nullstelle in $[g_n, g_{n+1}]$ und eine in $[g_{n+1}, g_{n+2}]$.

---

## Klassen

### `RiemannSiegelExtended`

| Methode | Beschreibung |
|---------|-------------|
| `theta(t)` | Siegelsche θ-Funktion |
| `Z(t, terms)` | Hardy-Z-Funktion via Riemann-Siegel-Formel |
| `error_bound(t)` | Obere Fehlerschranke $0.053 \cdot t^{-3/4} + 0.14 \cdot t^{-5/4}$ |
| `find_gram_point(n)` | Newton-Verfahren: $g_n$ mit $\theta(g_n) = n\pi$ |
| `gram_sign(n)` | Erwartetes Vorzeichen $(-1)^n$ |
| `check_gram_law(n)` | Prüft ob Gram-Gesetz bei $g_n$ gilt |
| `gram_statistics(n_max)` | Statistik: Erfolgsrate des Gram-Gesetzes |
| `find_zeros_in_interval(a,b)` | Nullstellen von Z(t) via Bisektionsverfahren |

### `GramBlocks`

| Methode | Beschreibung |
|---------|-------------|
| `identify_gram_blocks(n_start, n_end)` | Identifiziert Gram-Blöcke |
| `rosser_rule_holds(n)` | Prüft Rosser-Regel für $g_n$ |
| `gram_block_statistics(n_max)` | Statistik der Gram-Block-Breiten |

---

## Bekannte Nullstellen der ζ-Funktion

Die ersten nichttrivialen Nullstellen $\rho = 1/2 + it_k$:

| $k$ | $t_k$ (bekannt) |
|-----|-----------------|
| 1 | 14.134725141734694... |
| 2 | 21.022039638771554... |
| 3 | 25.010857580145688... |
| 4 | 30.424876125859513... |
| 5 | 32.935061587739189... |

---

## Hinweise zur Riemann-Hypothese

Die **Riemann-Hypothese** (1859, Conjecture — OFFEN) besagt:
> Alle nichttrivialen Nullstellen von $\zeta(s)$ liegen auf der kritischen Geraden $\text{Re}(s) = 1/2$.

In Termen der Z-Funktion: **Alle reellen Nullstellen von $Z(t)$ für $t > 0$** entsprechen
nichttrivialen Nullstellen auf der kritischen Geraden.

Numerisch verifiziert bis $t \approx 10^{13}$ (über $10^{13}$ Nullstellen).
