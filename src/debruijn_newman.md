# debruijn_newman.py – de Bruijn-Newman-Konstante Λ

**Autor**: Michael Fuhrmann
**Version**: 1.0
**Letzte Änderung**: 2026-03-12

---

## Mathematischer Hintergrund

### Definition der de Bruijn-Newman-Konstante

Die **de Bruijn-Newman-Konstante** $\Lambda$ ist definiert durch eine parametrische Familie von Funktionen:

$$H_t(x) = \int_0^\infty e^{tu^2} \cdot \Phi(u) \cdot \cos(xu) \, du$$

wobei die Kernfunktion:

$$\Phi(u) = \sum_{n=1}^\infty \left(2\pi^2 n^4 e^{9u} - 3\pi n^2 e^{5u}\right) \exp\left(-\pi n^2 e^{4u}\right)$$

**Definition von $\Lambda$**:
$$\Lambda = \inf\{ t \in \mathbb{R} : H_t \text{ hat nur reelle Nullstellen} \}$$

### Verbindung zur Riemann-Hypothese

**Spezialfall $t = 0$**: Es gilt
$$H_0(x) = \tfrac{1}{2} \xi\!\left(\tfrac{1}{2} + ix\right)$$

wo $\xi(s) = \tfrac{1}{2} s(s-1) \pi^{-s/2} \Gamma(s/2) \zeta(s)$ die **vollständige Riemann-Zeta-Funktion** ist.

**Die Äquivalenz** (de Bruijn 1950 / Newman 1976):
$$\boxed{\text{RH} \iff \Lambda \leq 0}$$

Da Rodgers-Tao 2018 bewiesen haben, dass $\Lambda \geq 0$, gilt:
$$\text{RH wahr} \iff \Lambda = 0$$

### Nullstellen von $H_0$

Die Nullstellen von $H_0$ entsprechen den imaginären Teilen der Riemann-Zeta-Nullstellen auf der kritischen Geraden $\text{Re}(s) = \frac{1}{2}$:

| $k$ | $\gamma_k$ (bekannt) |
|-----|---------------------|
| 1 | 14.134725141734... |
| 2 | 21.022039638771... |
| 3 | 25.010857580145... |
| 4 | 30.424876125859... |
| 5 | 32.935061587739... |

### de Bruijn-Satz (1950)

Wenn $H_{t_0}$ nur reelle Nullstellen hat, dann hat $H_t$ nur reelle Nullstellen für alle $t \geq t_0$. Dies macht $\Lambda$ zu einem **wohldefinierten** unteren Schranken-Wert.

**Beweisskizze**: $H_t(x) = H_0 \ast G_t$ (Faltung mit Gauß-Kern $G_t(x) \propto e^{-x^2/(4t)}$ für $t > 0$). Die Faltung mit einem Gauß-Kern verschiebt komplexe Nullstellen in Richtung der reellen Achse.

### Newman-Vermutung und Rodgers-Tao-Beweis

**Newman (1976)** vermutete: $\Lambda \geq 0$ (die RH ist "scharf").

**Rodgers & Tao (2018)** bewiesen $\Lambda \geq 0$ mit wahrscheinlichkeitstheoretischen Methoden, insbesondere durch Analyse von Korrelationen der GUE-Nullstellenverteilung.

**Konsequenz**: Falls RH wahr ist, dann $\Lambda = 0$ (nicht $\Lambda < 0$).

---

## Historische Schranken für $\Lambda$

| Jahr | Schranke | Quelle |
|------|----------|--------|
| 1950 | $\Lambda \leq \frac{1}{2}$ | de Bruijn |
| 1976 | Vermutung $\Lambda \geq 0$ | Newman |
| 2009 | $\Lambda \geq -2.7 \times 10^{-9}$ | Saouter, Gourdon, Demichel |
| 2018 | $\Lambda \geq 0$ (**bewiesen**) | Rodgers & Tao |
| 2020 | $0 \leq \Lambda \leq 0.22$ | Polymath15 |
| 2021 | $0 \leq \Lambda \leq 0.2$ | Platt & Trudgian |

---

## Klassen-API

### `DeBruijnNewman(praezision=50)`

Initialisiert mit mpmath-Präzision in Dezimalstellen.

### `phi_funktion(u, N_max=20) -> float`

Berechnet $\Phi(u)$ mit $N_{\max}$ Summanden. Für $u \geq 0$ dominiert $n=1$; bereits $N_{\max}=5$ reicht für 10+ Dezimalstellen.

### `H_t(x, t, N_max=30) -> float`

Berechnet $H_t(x)$ via adaptiver Gauss-Quadratur (mpmath) oder Simpson-Fallback. Integrationsbereich $[0, 5]$ (Φ ist danach $\approx 0$).

### `suche_nullstellen_H_t(t, x_bereich, N_max, schritte) -> List[float]`

Sucht Nullstellen via Vorzeichenwechsel + Bisektionsverfahren (50 Iterationen → ~15 Dezimalstellen).

### `ist_alle_reell(t, grenze, N_max, toleranz) -> bool`

Heuristische Prüfung: Hat $H_t$ im Bereich $[0, \text{grenze}]$ nur reelle Nullstellen?

### `schranke_lambda_oben(kandidat_t) -> Dict`

Gibt numerische Evidenz für $\Lambda \leq \text{kandidat\_t}$. Hinweis: Kein rigoroser Beweis.

### `riemann_zusammenhang() -> Dict`

Dokumentiert $\text{RH} \iff \Lambda \leq 0$ formal.

### `newman_vermutung_status() -> Dict`

Status-Bericht: $\Lambda \geq 0$ bewiesen (Rodgers-Tao 2018).

### `bekannte_schranken_history() -> List[Dict]`

Chronologische Liste aller bekannten Schranken für $\Lambda$.

### `nullstellen_H0_berechnen(x_max, N_max) -> List[float]`

Berechnet Nullstellen von $H_0(x)$ (= Riemann-Zeta-Nullstellen).

### `vergleiche_mit_riemann_nullstellen(nullstellen) -> Dict`

Vergleicht numerisch gefundene Nullstellen mit bekannten $\gamma_k$-Werten.

### `visualisiere_H0(x_max, N_max, dateiname) -> bool`

Plottet $H_0(x)$ mit markierten Riemann-Nullstellen (benötigt matplotlib).

---

## Numerische Einschränkungen

- **Integrationsgenauigkeit**: Die Quadratur von $H_t(x)$ für großes $x$ ist durch Oszillation des $\cos(xu)$-Terms begrenzt.
- **N_max-Konvergenz**: Für $u \ll 1$ ist $\Phi(u)$ hauptsächlich durch $n=1$ bestimmt; für $u \gg 1$ fallen alle Terme sehr schnell ab.
- **Riemann-Nullstellen**: Die numerischen Nullstellen von $H_0$ können von den exakten Werten abweichen (Fehler ~0.1 bis 1.0 je nach N_max und Schrittweite).

---

## Quellen

- N.G. de Bruijn (1950): *The roots of trigonometric integrals*, Duke Math. J.
- C.M. Newman (1976): *Fourier transforms with only real zeros*, Proc. AMS
- B. Rodgers, T. Tao (2018): *The de Bruijn-Newman constant is non-negative*, arXiv:1801.05914
- D. Platt, T. Trudgian (2021): *The Riemann hypothesis is true up to $3 \times 10^{12}$*
- Polymath15 project (2020): *Effective approximation of heat flow...*
