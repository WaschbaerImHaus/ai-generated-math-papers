# pde.py — Partielle Differentialgleichungen

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-10
**Abhängigkeiten:** `numpy`, `scipy.sparse`, `scipy.sparse.linalg`, `scipy.linalg`

---

## Überblick

Das Modul `pde.py` implementiert numerische und analytische Lösungsverfahren für klassische Partielle Differentialgleichungen (PDEs). Es deckt die wichtigsten Gleichungstypen und Diskretisierungsmethoden ab, die in Physik, Ingenieurwesen und Mathematik Anwendung finden.

---

## 1. Klassifikation von PDEs 2. Ordnung

Jede lineare PDE 2. Ordnung in zwei Variablen lautet:

$$A u_{xx} + B u_{xy} + C u_{yy} + \text{(Terme niedriger Ordnung)} = 0$$

Die Diskriminante $D = B^2 - 4AC$ bestimmt den Typ:

| Bedingung | Typ | Beispiel |
|-----------|-----|---------|
| $D < 0$ | **Elliptisch** | Laplace: $\Delta u = 0$ |
| $D = 0$ | **Parabolisch** | Wärme: $u_t = \alpha^2 u_{xx}$ |
| $D > 0$ | **Hyperbolisch** | Welle: $u_{tt} = c^2 u_{xx}$ |

---

## 2. Wärmeleitungsgleichung

$$u_t = \alpha^2 u_{xx}, \quad x \in [0,L],\; t > 0$$

mit Dirichlet-Randwerten $u(0,t) = u(L,t) = 0$ und Anfangsbedingung $u(x,0) = u_0(x)$.

### 2.1 Explizites FTCS-Schema (Forward-Time Central-Space)

$$u_i^{n+1} = u_i^n + r\bigl(u_{i+1}^n - 2u_i^n + u_{i-1}^n\bigr), \quad r = \frac{\alpha^2 \Delta t}{\Delta x^2}$$

**Stabilitätsbedingung (von Neumann):**
$$r = \frac{\alpha^2 \Delta t}{\Delta x^2} \leq \frac{1}{2}$$

- Einfach zu implementieren, aber stark zeitschrittbeschränkt.

### 2.2 Implizites BTCS-Schema (Backward-Time Central-Space)

$$(I + r K)\, u^{n+1} = u^n$$

wobei $K$ die diskrete zweite Ableitung (Tridiagonalmatrix) ist.

- **Unbedingt stabil** für alle $r > 0$.
- Erfordert Lösung eines tridiagonalen LGS in jedem Schritt ($O(n)$).

### 2.3 Crank-Nicolson-Schema ($\theta = \tfrac{1}{2}$)

$$\Bigl(I + \tfrac{r}{2} K\Bigr)\, u^{n+1} = \Bigl(I - \tfrac{r}{2} K\Bigr)\, u^n$$

- **Zweite Ordnung** in Zeit und Raum: $O(\Delta t^2 + \Delta x^2)$.
- **Unbedingt stabil**.
- Bevorzugtes Verfahren für Wärmeleitungsgleichungen.

### 2.4 Analytische Fourier-Reihen-Lösung

Für $u(x,0) = u_0(x)$, $u(0,t) = u(L,t) = 0$:

$$u(x,t) = \sum_{n=1}^{\infty} b_n \sin\!\Bigl(\frac{n\pi x}{L}\Bigr) e^{-(n\pi/L)^2 t}$$

mit Fourier-Koeffizienten:
$$b_n = \frac{2}{L} \int_0^L u_0(x) \sin\!\Bigl(\frac{n\pi x}{L}\Bigr)\,dx$$

---

## 3. Wellengleichung

$$u_{tt} = c^2 u_{xx}, \quad x \in [0,L],\; t > 0$$

mit $u(0,t) = u(L,t) = 0$, $u(x,0) = u_0(x)$, $u_t(x,0) = v_0(x)$.

### 3.1 Leapfrog-Schema

$$u_i^{n+1} = 2u_i^n - u_i^{n-1} + \lambda^2\bigl(u_{i+1}^n - 2u_i^n + u_{i-1}^n\bigr)$$

**CFL-Bedingung (Courant-Friedrichs-Lewy):**
$$\lambda = \frac{c\,\Delta t}{\Delta x} \leq 1$$

- Bei CFL-Verletzung: **instabiles** Verhalten.

### 3.2 Analytische Lösung (Fourier)

$$u(x,t) = \sum_{n=1}^{\infty} \Bigl[a_n \cos\!\Bigl(\frac{n\pi c t}{L}\Bigr) + b_n \sin\!\Bigl(\frac{n\pi c t}{L}\Bigr)\Bigr] \sin\!\Bigl(\frac{n\pi x}{L}\Bigr)$$

mit:
$$a_n = \frac{2}{L}\int_0^L u_0(x)\sin\!\Bigl(\frac{n\pi x}{L}\Bigr)dx, \quad b_n = \frac{2}{n\pi c}\int_0^L v_0(x)\sin\!\Bigl(\frac{n\pi x}{L}\Bigr)dx$$

### 3.3 Dispersionsrelation

Die Wellengleichung ist **nicht-dispersiv**:
$$\omega = ck$$

Phasengeschwindigkeit = Gruppengeschwindigkeit = $c$.

---

## 4. Laplace- und Poisson-Gleichung

### 4.1 Laplace-Gleichung

$$\Delta u = u_{xx} + u_{yy} = 0 \quad \text{auf } [0,1]\times[0,1]$$

mit Dirichlet-Randbedingungen.

**Numerik:** Successive Over-Relaxation (SOR) mit optimalem Parameter:
$$\omega_{\rm opt} = \frac{2}{1 + \sin(\pi/N)}$$

Iterationsschritt:
$$u_{i,j}^{\rm new} = (1-\omega)\,u_{i,j} + \frac{\omega}{4}\bigl(u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}\bigr)$$

### 4.2 Poisson-Gleichung

$$-\Delta u = f(x,y) \quad \text{auf } [0,1]\times[0,1], \quad u = 0 \text{ am Rand}$$

**Numerik:** Direkte Lösung via sparse LGS (5-Punkt-Stern-Diskretisierung):

$$-u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1} + 4u_{i,j} = \Delta x^2 f_{i,j}$$

### 4.3 Greens Funktion von Laplace in 2D

Die Fundamentallösung:
$$G(x,y;\,x_0,y_0) = -\frac{1}{2\pi} \ln\!\bigl(\|(x,y)-(x_0,y_0)\|\bigr)$$

erfüllt $-\Delta_x G = \delta(x-x_0)\delta(y-y_0)$.

### 4.4 Eigenschaften harmonischer Funktionen

Harmonische Funktionen ($\Delta u = 0$) erfüllen:
- **Mittelwerteigenschaft:** $u(x_0) = \frac{1}{|\partial B_r|}\int_{\partial B_r} u\,dS$
- **Maximum-Prinzip:** Extremwerte werden nur am Rand angenommen

---

## 5. Schrödinger-Gleichung

### 5.1 Stationäre Schrödinger-Gleichung

$$-\frac{\hbar^2}{2m}\psi'' + V(x)\psi = E\psi \quad (\hbar = m = 1)$$

**Numerik:** Tridiagonales Eigenwertproblem (LAPACK `eigh_tridiagonal`):

$$-\psi_{i+1} + 2\psi_i - \psi_{i-1} = \Delta x^2\,(E - V_i)\,\psi_i$$

### 5.2 Harmonischer Oszillator

Potential: $V(x) = \tfrac{1}{2}x^2$

Analytische Lösung:
$$E_n = n + \tfrac{1}{2}, \quad n = 0, 1, 2, \ldots$$
$$\psi_n(x) = (2^n n! \sqrt{\pi})^{-1/2}\, H_n(x)\, e^{-x^2/2}$$

wobei $H_n$ die physikalischen Hermite-Polynome sind.

### 5.3 Zeitabhängige Schrödinger-Gleichung

$$i\hbar \frac{\partial\psi}{\partial t} = \hat{H}\psi = \Bigl[-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x)\Bigr]\psi$$

**Crank-Nicolson-Schema** (unitär, norm-erhaltend):

$$\Bigl(I + \tfrac{i\Delta t}{2}H\Bigr)\psi^{n+1} = \Bigl(I - \tfrac{i\Delta t}{2}H\Bigr)\psi^n$$

Normerhaltung: $\|\psi^{n+1}\|^2 = \|\psi^n\|^2$.

---

## 6. Finite-Elemente-Methode (FEM) 1D

Für $-u'' = f(x)$, $u(a) = u(b) = 0$:

**Schwache Formulierung (Galerkin):**
$$\int_a^b u'v'\,dx = \int_a^b f\,v\,dx \quad \forall v \in H_0^1(a,b)$$

**Hutfunktionen** $\phi_k(x)$: stückweise linear, $\phi_k(x_j) = \delta_{kj}$.

**Steifigkeitsmatrix** (Element $[x_k, x_{k+1}]$, $h = \Delta x$):
$$K_{\rm loc} = \frac{1}{h}\begin{pmatrix}1 & -1\\-1 & 1\end{pmatrix}$$

**Lastvektor:**
$$F_{\rm loc} = \frac{h}{6}\begin{pmatrix}2f(x_k)+f(x_{k+1})\\f(x_k)+2f(x_{k+1})\end{pmatrix}$$

**Konvergenz:** $\|u_h - u\|_{L^2} = O(h^2)$ für glatte Lösungen.

---

## 7. Methode der Charakteristiken

Für $a\,u_x + b\,u_t = c\,u$:

Charakteristiken: $\frac{dx}{ds} = a$, $\frac{dt}{ds} = b$ → Geraden $x = x_0 + \frac{a}{b}t$.

Entlang einer Charakteristik:
$$u(x,t) = u_0\!\Bigl(x - \frac{a}{b}t\Bigr) \cdot e^{(c/b)\,t}$$

---

## 8. Viskose Burgers-Gleichung

$$u_t + u\,u_x = \nu\,u_{xx}$$

Modellproblem für:
- Nichtlineare Advektion (Schockbildung für $\nu \to 0$)
- Viskose Dissipation

**Upwind-Schema (explizit, periodisch):**
$$u_i^{n+1} = u_i^n + \Delta t\Bigl(-u_i\,[u_x]_i^{\rm upwind} + \nu\,\frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2}\Bigr)$$

---

## Numerische Stabilität (Zusammenfassung)

| Methode | Schema | Stabilitätsbedingung | Ordnung |
|---------|--------|----------------------|---------|
| Wärme explizit | FTCS | $r = \alpha^2\Delta t/\Delta x^2 \leq 1/2$ | $O(\Delta t + \Delta x^2)$ |
| Wärme implizit | BTCS | unbedingt stabil | $O(\Delta t + \Delta x^2)$ |
| Crank-Nicolson | $\theta=1/2$ | unbedingt stabil | $O(\Delta t^2 + \Delta x^2)$ |
| Welle Leapfrog | - | CFL: $c\Delta t/\Delta x \leq 1$ | $O(\Delta t^2 + \Delta x^2)$ |
| Schrödinger CN | unitär | unbedingt stabil | $O(\Delta t^2 + \Delta x^2)$ |
| FEM 1D | Galerkin | - | $O(h^2)$ in $L^2$ |
