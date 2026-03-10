# tensor_geometry.py — Tensor- und Differentialgeometrie

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Modul:** `src/tensor_geometry.py`

---

## Übersicht

Das Modul `tensor_geometry.py` implementiert die mathematischen Grundlagen der **Tensor- und Differentialgeometrie** — jenes Teilgebiet der Mathematik, das die Geometrie allgemeiner gekrümmter Räume beschreibt. Es bildet das mathematische Fundament der **Allgemeinen Relativitätstheorie** (Einstein 1915) und ist ein aktives Forschungsgebiet der modernen Mathematik und Physik.

---

## 1. Tensoren

### 1.1 Was ist ein Tensor?

Ein Tensor ist eine **multilineare Abbildung**, die Vektoren und Kovektoren auf reelle Zahlen abbildet, und die sich bei Koordinatenwechsel nach bestimmten Transformationsregeln verhält.

Klassifikation nach **Rang $(p,q)$**:

| Rang | Bezeichnung | Beispiel |
|------|-------------|---------|
| $(0,0)$ | Skalar | Temperatur, Druck |
| $(1,0)$ | Kontravariantervektor | Geschwindigkeit $v^i$ |
| $(0,1)$ | Kovarianter Vektor / 1-Form | Gradient $\nabla f$ |
| $(1,1)$ | Lineare Abbildung | Spannungstensor |
| $(0,2)$ | Bilinearform | Metriktensor $g_{ij}$ |
| $(2,0)$ | Bikontravariant | Inverse Metrik $g^{ij}$ |
| $(1,3)$ | Rang-4-Tensor | Riemann-Tensor $R^i{}_{jkl}$ |

### 1.2 Einstein-Summenkonvention

Über doppelt auftretende Indizes (einmal oben, einmal unten) wird automatisch summiert:

$$T^{ij} S_{jk} = \sum_{j=0}^{n-1} T^{ij} S_{jk}$$

### 1.3 Koordinatentransformation

Bei einem Koordinatenwechsel $x^i \mapsto x^{i'}$ transformiert sich ein $(1,1)$-Tensor als:

$$T^{i'}{}_{j'} = \frac{\partial x^{i'}}{\partial x^i} \cdot \frac{\partial x^j}{\partial x^{j'}} \cdot T^i{}_j$$

Kontravariante Indizes transformieren mit $\partial x'/\partial x$, kovariante mit $\partial x/\partial x'$.

### 1.4 Tensoroperationen

#### Äußeres Produkt (Tensorprodukt) $\otimes$

$$
(T \otimes S)^{i_1 i_2}{}_{j_1 j_2} = T^{i_1}{}_{j_1} \cdot S^{i_2}{}_{j_2}
$$

Rang: $(p_1+p_2,\; q_1+q_2)$.

#### Kontraktion (Verjüngung)

Setzt oberen Index $i$ gleich unterem Index $j$ und summiert:

$$
C(T)^{i_1\ldots\hat{i}\ldots}{}_{j_1\ldots\hat{j}\ldots} = \sum_k T^{i_1\ldots k\ldots}{}_{j_1\ldots k\ldots}
$$

Rang sinkt um $(1,1)$. Beispiel: Spur einer Matrix $= T^i{}_i$.

#### Symmetrisierung und Antisymmetrisierung

$$
T_{(ij)} = \frac{1}{2}(T_{ij} + T_{ji}), \qquad T_{[ij]} = \frac{1}{2}(T_{ij} - T_{ji})
$$

Zerlegung: $T_{ij} = T_{(ij)} + T_{[ij]}$ (eindeutig).

---

## 2. Metriktensor

### 2.1 Definition

Der Metriktensor $g_{ij}$ ist ein $(0,2)$-Tensor, der an jedem Punkt einer Mannigfaltigkeit ein inneres Produkt auf dem Tangentialraum definiert.

**Linienelement:**
$$ds^2 = g_{ij}\, dx^i\, dx^j$$

**Eigenschaften:**
- Symmetrisch: $g_{ij} = g_{ji}$
- Nicht-degeneriert: $\det(g) \neq 0$
- Riemannsch: positiv definit (alle Eigenwerte $> 0$)
- Lorentzsich (ART): Signatur $(-,+,+,+)$

### 2.2 Indexhebung und -senkung

Mit der inversen Metrik $g^{ij}$ (sodass $g^{ik} g_{kj} = \delta^i_j$):

$$T^i = g^{ij} T_j \quad \text{(Indexhebung, musical isomorphism } \sharp\text{)}$$
$$T_i = g_{ij} T^j \quad \text{(Indexsenkung, musical isomorphism } \flat\text{)}$$

### 2.3 Volumenelement

$$dV = \sqrt{|\det(g)|}\; dx^1 \wedge dx^2 \wedge \cdots \wedge dx^n$$

---

## 3. Christoffel-Symbole

### 3.1 Definition

Die **Christoffel-Symbole** $\Gamma^k{}_{ij}$ (auch: Levi-Civita-Verbindungskoeffizienten) beschreiben, wie Koordinatenbasisvektoren sich von Punkt zu Punkt ändern:

$$
\Gamma^k{}_{ij} = \frac{1}{2}\, g^{kl}\!\left(\partial_i g_{jl} + \partial_j g_{il} - \partial_l g_{ij}\right)
$$

### 3.2 Eigenschaften

- **Kein Tensor** (transformiert sich anders als Tensoren)
- **Symmetrisch** in unteren Indizes: $\Gamma^k{}_{ij} = \Gamma^k{}_{ji}$ (torsionsfrei)
- **Flacher Raum** (kart. Koordinaten): $\Gamma^k{}_{ij} = 0$ überall
- Anzahl: $n^2(n+1)/2$ unabhängige Komponenten in $n$ Dimensionen

### 3.3 Bekannte Werte (Sphäre, $r=1$)

Sphärische Koordinaten $(\theta, \phi)$:

$$\Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta, \qquad \Gamma^\phi{}_{\theta\phi} = \frac{\cos\theta}{\sin\theta}$$

Alle anderen Christoffel-Symbole der Sphäre verschwinden.

---

## 4. Riemannscher Krümmungstensor

### 4.1 Definition

$$R^l{}_{kij} = \partial_i \Gamma^l{}_{jk} - \partial_j \Gamma^l{}_{ik} + \Gamma^l{}_{im}\Gamma^m{}_{jk} - \Gamma^l{}_{jm}\Gamma^m{}_{ik}$$

### 4.2 Geometrische Bedeutung

**Anschauung:** Transportiert man einen Vektor parallel entlang einer kleinen, geschlossenen Kurve auf einer gekrümmten Fläche, kehrt er **gedreht** zurück. Der Drehwinkel ist proportional zur Fläche der Kurve und zum Riemann-Tensor.

Im flachen Raum: $R^l{}_{kij} = 0$ überall (Holonomie trivial).

### 4.3 Symmetrien

$$R_{lkij} = -R_{lkji}, \qquad R_{lkij} = -R_{klij}, \qquad R_{lkij} = R_{ijlk}$$

**Bianchi-Identität:** $R^l{}_{k[ij;m]} = 0$ (kovariante zyklische Summe = 0).

Unabhängige Komponenten: $\frac{n^2(n^2-1)}{12}$ (in 4D: 20).

---

## 5. Ricci-Tensor und Ricci-Skalar

### 5.1 Ricci-Tensor

Kontraktion des Riemann-Tensors:

$$R_{ij} = R^k{}_{ikj}$$

- Symmetrisch: $R_{ij} = R_{ji}$
- Misst Abweichung des Volumens einer geodätischen Kugel vom flachen Fall

### 5.2 Ricci-Skalar

$$R = g^{ij} R_{ij}$$

Vollständige Kontraktion — eine einzige reelle Zahl.

### 5.3 Gaußsche Krümmung (2D)

$$K = \frac{R}{2} \quad \text{(nur für 2D-Flächen)}$$

---

## 6. Klassische Mannigfaltigkeiten — Krümmungsübersicht

| Mannigfaltigkeit | Koordinaten | Metrik $ds^2$ | Gaußsche Krümmung $K$ |
|-----------------|-------------|--------------|----------------------|
| **Ebene** $\mathbb{R}^2$ | $(x,y)$ | $dx^2 + dy^2$ | $0$ |
| **Sphäre** $S^2_r$ | $(\theta,\phi)$ | $r^2(d\theta^2 + \sin^2\!\theta\, d\phi^2)$ | $1/r^2$ |
| **Hyperb. Ebene** $H^2$ | $(x,y), y>0$ | $(dx^2+dy^2)/y^2$ | $-1$ |
| **Torus** $T^2$ | $(\theta,\phi)$ | $r^2 d\theta^2 + (R+r\cos\theta)^2 d\phi^2$ | $\frac{\cos\theta}{r(R+r\cos\theta)}$ |
| **Sattelfläche** $z=xy$ | $(x,y)$ | $(1+y^2)dx^2 + 2xy\,dxdy + (1+x^2)dy^2$ | $-\frac{1}{(1+x^2+y^2)^2}$ |
| **Schwarzschild** | $(t,r)$ | $-(1-r_s/r)c^2dt^2 + \frac{dr^2}{1-r_s/r}$ | (4D, kein $K$) |

---

## 7. Geodäten

### 7.1 Geodätengleichung

Eine **Geodäte** ist die Verallgemeinerung einer „geraden Linie" auf einer gekrümmten Mannigfaltigkeit — die kürzeste (oder längste) Verbindung zweier Punkte.

$$\frac{d^2 x^k}{dt^2} + \Gamma^k{}_{ij}\, \frac{dx^i}{dt}\, \frac{dx^j}{dt} = 0$$

Als System erster Ordnung (mit $v^k = \dot{x}^k$):

$$\dot{x}^k = v^k, \qquad \dot{v}^k = -\Gamma^k{}_{ij}\, v^i v^j$$

**Implementierung:** Runge-Kutta 4. Ordnung.

### 7.2 Bogenlänge

$$L = \int_0^T \sqrt{g_{ij}(x(t))\, \dot{x}^i \dot{x}^j}\; dt$$

### 7.3 Beispiele für Geodäten

| Mannigfaltigkeit | Geodäte |
|-----------------|---------|
| Ebene | Gerade Linie |
| Sphäre | Großkreis |
| Hyperb. Ebene | Halbkreise und vertikale Geraden |
| Torus | Linienwicklungen |

---

## 8. Differentialformen

### 8.1 Äußeres Produkt (Wedge-Produkt)

Für 1-Formen $\alpha, \beta$:

$$(\alpha \wedge \beta)_{ij} = \alpha_i \beta_j - \alpha_j \beta_i$$

Eigenschaften:
- Antisymmetrisch: $\alpha \wedge \beta = -\beta \wedge \alpha$
- $\alpha \wedge \alpha = 0$
- Bilinear und assoziativ

### 8.2 Äußere Ableitung

Für ein Skalarfeld $f$ (0-Form):
$$df = \frac{\partial f}{\partial x^i}\, dx^i \quad \text{(Gradient als 1-Form)}$$

Für eine 1-Form $\alpha$:
$$(d\alpha)_{ij} = \partial_i \alpha_j - \partial_j \alpha_i$$

Eigenschaft: $d^2 = 0$ (zweimalige Ableitung ist null).

### 8.3 Hodge-Stern-Operator

$$* : \Omega^k(M) \to \Omega^{n-k}(M)$$

Für 2D mit Einheitsmetrik: $*(dx^1) = dx^2$, $*(dx^2) = -dx^1$.

Eigenschaft: $**\omega = (-1)^{k(n-k)}\omega$ (für Riemannsche Metriken).

### 8.4 Lie-Ableitung

$$\mathcal{L}_V f = V^i \partial_i f \quad \text{(direktionale Ableitung in Richtung } V\text{)}$$

Rate der Änderung von $f$ entlang des Flusses von $V$.

---

## 9. Einstein-Gleichungen

### 9.1 Einstein-Tensor

$$G_{ij} = R_{ij} - \frac{1}{2} g_{ij} R$$

Eigenschaften:
- Symmetrisch: $G_{ij} = G_{ji}$
- Divergenzfrei: $\nabla^i G_{ij} = 0$ (Bianchi-Identität → Energieerhaltung!)

### 9.2 Feldgleichungen der ART

$$G_{ij} = \frac{8\pi G}{c^4} T_{ij}$$

- $G_{ij}$: Einstein-Tensor (Geometrie der Raumzeit)
- $T_{ij}$: Energie-Impuls-Tensor (Materie und Energie)
- $G$: Gravitationskonstante, $c$: Lichtgeschwindigkeit

**Vakuumgleichungen** ($T_{ij} = 0$): $G_{ij} = 0 \Leftrightarrow R_{ij} = 0$ (Ricci-flach).

### 9.3 Schwarzschild-Metrik

Exakte Lösung der Vakuumgleichungen für eine kugelsymmetrische Masse $M$:

$$ds^2 = -\!\left(1-\frac{r_s}{r}\right)c^2\, dt^2 + \frac{dr^2}{1-r_s/r} + r^2\, d\Omega^2$$

**Schwarzschild-Radius:**
$$r_s = \frac{2GM}{c^2}$$

Für die Sonne: $r_s \approx 3\,\text{km}$ (die Sonne ist viel größer → kein schwarzes Loch).

---

## 10. Hilfstensoren

### 10.1 Levi-Civita-Symbol

$$\varepsilon_{i_1 \ldots i_n} = \begin{cases} +1 & \text{gerade Permutation von } (1,\ldots,n) \\ -1 & \text{ungerade Permutation} \\ 0 & \text{sonst (wiederholter Index)} \end{cases}$$

Verwendung: Determinante, Kreuzprodukt, Volumenform.

### 10.2 Kronecker-Delta

$$\delta^i{}_j = \begin{cases} 1 & i = j \\ 0 & i \neq j \end{cases}$$

### 10.3 Paralleltransport

Transportiert einen Vektor $v^k$ entlang einer Kurve $x(t)$, ohne ihn zu „drehen":

$$\frac{Dv^k}{dt} = \frac{dv^k}{dt} + \Gamma^k{}_{ij}\, v^i\, \dot{x}^j = 0$$

Euler-Diskretisierung: $\Delta v^k = -\Gamma^k{}_{ij}\, v^i\, \Delta x^j$.

---

## 11. Verbindungen zu anderen Gebieten

### Zur Physik

| Physikalisches Konzept | Mathematische Struktur |
|-----------------------|----------------------|
| Raumzeit-Krümmung | Riemannscher Krümmungstensor $R^l{}_{kij}$ |
| Gravitationsfeld | Christoffel-Symbole $\Gamma^k{}_{ij}$ |
| Energie-Impuls | Tensor $T_{ij}$ ($(0,2)$-Tensor) |
| Geodäte (Lichtweg) | Lösung der Geodätengleichung |
| Elektromagnetismus | Faraday-Tensor $F_{ij}$ (antisymm. $(0,2)$-Tensor) |
| Quantenfeldtheorie | Spinoren, Clifford-Algebren (Erweiterung) |
| Yang-Mills-Theorie | Verbindungsformen auf Faserbündeln |

### Zur Mathematik

| Mathematisches Konzept | Verbindung |
|-----------------------|----------|
| Differentialtopologie | Mannigfaltigkeiten als Grundobjekte |
| Algebraische Topologie | Euler-Charakteristik via Gaußschen Bonnet |
| Komplexe Analysis | Kähler-Mannigfaltigkeiten (komplexe Struktur + Metrik) |
| Funktionalanalysis | Differentialoperatoren auf Tensorfeldern |
| Kategorientheorie | Natürliche Transformationen zwischen Tensorfeldern |

---

## 12. Gaußscher Bonnet-Satz

Verbindet globale Topologie mit lokaler Geometrie:

$$\int_M K\, dA = 2\pi\, \chi(M)$$

- $K$: Gaußsche Krümmung
- $dA = \sqrt{|\det(g)|}\, dx^1\, dx^2$: Flächenelement
- $\chi(M)$: Euler-Charakteristik ($\chi(S^2) = 2$, $\chi(T^2) = 0$)

Konsequenz: Man kann die Topologie einer Fläche allein durch ihre Krümmung bestimmen!

---

## 13. Implementierungsdetails

### Numerische Ableitung

Alle Ableitungen von Metriken und Christoffel-Symbolen werden numerisch berechnet via **zentralem Differenzenquotienten** (Genauigkeit $O(h^2)$):

$$\frac{\partial f}{\partial x^k} \approx \frac{f(x + h\,e_k) - f(x - h\,e_k)}{2h}$$

Typische Schrittweiten: $h = 10^{-5}$ (Christoffel), $h = 10^{-4}$ (Riemann-Tensor).

### Runge-Kutta 4 für Geodäten

Die Geodätengleichung wird als System erster Ordnung formuliert und mit dem klassischen RK4-Verfahren (4 Stufen, Ordnung 4) integriert. Lokaler Fehler: $O((\Delta t)^5)$.

### Singularitäten

- **Sphäre bei $\theta = 0, \pi$**: Koordinaten-Singularität (kein Fehler in der Geometrie)
- **Schwarzschild bei $r = r_s$**: Koordinaten-Singularität (Ereignishorizont); echte Singularität bei $r = 0$
- **Hyperbolische Ebene bei $y = 0$**: Rand der Halbebene; kleine Regularisierung $y \to \max(y, 10^{-12})$
