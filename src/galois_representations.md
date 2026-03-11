# Galois-Darstellungen — Mathematische Dokumentation

**Datei:** `src/galois_representations.py`
**Autor:** Michael Fuhrmann
**Version:** 1.0
**Stand:** 2026-03-11

---

## Überblick

Das Modul implementiert die Theorie der **ℓ-adischen Galois-Darstellungen**, eine der
zentralen Brücken zwischen algebraischer Zahlentheorie, arithmetischer Geometrie und
dem Langlands-Programm.

Eine **Galois-Darstellung** ist ein stetiger Gruppenhomomorphismus

$$\rho: \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q}) \longrightarrow \mathrm{GL}_n(\mathbb{Z}_\ell)$$

wobei $\mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})$ die absolute Galois-Gruppe von $\mathbb{Q}$
(mit profiniter Topologie) und $\mathbb{Z}_\ell$ die $\ell$-adischen ganzen Zahlen sind.

---

## 1. Abstrakte Basisklasse `GaloisRepresentation`

### Mathematischer Hintergrund

Die absolute Galois-Gruppe $G_\mathbb{Q} = \mathrm{Gal}(\bar{\mathbb{Q}}/\mathbb{Q})$ ist eine
profinite Gruppe — der projektive Limes aller endlichen Galois-Gruppen $\mathrm{Gal}(K/\mathbb{Q})$
über alle endlichen galoisschen Erweiterungen $K/\mathbb{Q}$.

Für fast alle Primzahlen $p$ (alle außer dem Führer der Darstellung) gibt es ein
**Frobenius-Element** $\mathrm{Frob}_p \in G_\mathbb{Q}$, das durch seine Wirkung auf
$p^n$-ten Einheitswurzeln charakterisiert ist:

$$\mathrm{Frob}_p(\zeta) = \zeta^p$$

Die **Frobenius-Spur** $\mathrm{tr}(\rho(\mathrm{Frob}_p))$ und das **charakteristische Polynom**
$\det(1 - \rho(\mathrm{Frob}_p) \cdot X)$ sind fundamentale arithmetische Invarianten.

### Lokaler L-Faktor

$$L_p(\rho, s) = \det\bigl(1 - \rho(\mathrm{Frob}_p) \cdot p^{-s}\bigr)^{-1}$$

Das Euler-Produkt $L(\rho, s) = \prod_p L_p(\rho, s)$ konvergiert für $\mathrm{Re}(s) \gg 1$.

---

## 2. Tate-Modul `TateModuleRepresentation`

### Konstruktion

Für eine elliptische Kurve $E/\mathbb{Q}$ und eine Primzahl $\ell$ ist der **ℓ-adische Tate-Modul**

$$T_\ell(E) = \varprojlim_n E[\ell^n]$$

der projektive Limes der $\ell^n$-Torsionspunkte bezüglich der Multiplikations-mit-$\ell$-Abbildung.

Als abelsche Gruppe gilt $T_\ell(E) \cong \mathbb{Z}_\ell^2$ (für $\ell \neq \mathrm{char}$).

### Galois-Wirkung

Die kanonische Wirkung von $G_\mathbb{Q}$ auf Torsionspunkten induziert eine Darstellung

$$\rho_{E,\ell}: G_\mathbb{Q} \longrightarrow \mathrm{GL}_2(\mathbb{Z}_\ell)$$

### Frobenius und Hasse-Weil

Bei **guter Reduktion** bei $p$ gilt für das charakteristische Polynom:

$$\det\bigl(1 - \mathrm{Frob}_p \cdot X \mid T_\ell(E)\bigr) = 1 - a_p X + p X^2$$

wobei $a_p = p + 1 - \#E(\mathbb{F}_p)$ die **Hasse-Weil-Spurformel** ist.

### Hasse-Schranke

$$|a_p| \leq 2\sqrt{p}$$

(bewiesen durch Hasse 1933 für elliptische Kurven, allgemein durch Deligne als "Weil II").

Äquivalent: Die Frobenius-Eigenwerte $\alpha, \beta$ mit $\alpha + \beta = a_p$ und $\alpha \beta = p$
erfüllen $|\alpha| = |\beta| = \sqrt{p}$ (liegen auf dem Kreis mit Radius $\sqrt{p}$).

### Beispiel: $E: y^2 = x^3 - x$ (Führer $N = 32$)

Diese Kurve hat besondere Symmetrie ($j$-Invariante $= 1728$, CM durch $\mathbb{Z}[i]$):
- Für $p \equiv 3 \pmod{4}$: $a_p = 0$ (da $p$ inert in $\mathbb{Z}[i]$)
- Für $p \equiv 1 \pmod{4}$: $a_p = 2a$ wo $p = a^2 + b^2$, $a \equiv 1 \pmod{2}$

---

## 3. Zyklotomischer Charakter `CyclotomicCharacter`

### Definition

Der **zyklotomische Charakter** $\chi_\ell: G_\mathbb{Q} \to \mathbb{Z}_\ell^*$ ist definiert durch
die Wirkung auf $\ell$-Potenz-Einheitswurzeln: für $\sigma \in G_\mathbb{Q}$ und primitive $\ell^n$-te
Einheitswurzel $\zeta$ gilt

$$\sigma(\zeta) = \zeta^{\chi_\ell(\sigma)}$$

### Frobenius-Auswertung

$$\chi_\ell(\mathrm{Frob}_p) = p \in \mathbb{Z}_\ell^*$$

In der Implementierung: $\chi_\ell(\mathrm{Frob}_p) = p \bmod \ell$ (Reduktion modulo $\ell$).

### Parität

$\chi_\ell$ ist **ungerade** (odd): Die komplexe Konjugation $c$ wirkt durch $\zeta \mapsto \zeta^{-1}$,
also $\chi_\ell(c) = -1 \in \mathbb{Z}_\ell^*$.

### Tensorpotenzen

$\chi_\ell^n(\mathrm{Frob}_p) = p^n \bmod \ell$

---

## 4. Dirichlet-Charakter-Darstellung `DirichletCharacterRepresentation`

### Kronecker-Weber und Klassenkörpertheorie

Nach dem **Kronecker-Weber-Theorem** ist jede abelsche Erweiterung $K/\mathbb{Q}$ in einem
zyklotomischen Körper $\mathbb{Q}(\zeta_n)$ enthalten.

Die zugehörigen 1-dimensionalen Galois-Darstellungen

$$\rho_\chi: \mathrm{Gal}(\mathbb{Q}(\zeta_q)/\mathbb{Q}) \xrightarrow{\;\sim\;} (\mathbb{Z}/q\mathbb{Z})^* \xrightarrow{\;\chi\;} \mathbb{C}^*$$

entsprechen genau den **Dirichlet-Charakteren** $\chi \pmod{q}$.

### Frobenius-Wirkung

$$\rho_\chi(\mathrm{Frob}_p) = \chi(p) \quad \text{für } p \nmid q$$

### L-Funktion

$$L(\rho_\chi, s) = L(s, \chi) = \sum_{n=1}^\infty \frac{\chi(n)}{n^s} = \prod_{p \nmid q} \frac{1}{1 - \chi(p) p^{-s}}$$

---

## 5. Symmetrische Potenz `SymmetricPowerRepresentation`

### Definition

Für eine 2-dimensionale Darstellung $\rho: G_\mathbb{Q} \to \mathrm{GL}_2(\mathbb{Z}_\ell)$ ist

$$\mathrm{Sym}^k(\rho): G_\mathbb{Q} \longrightarrow \mathrm{GL}_{k+1}(\mathbb{Z}_\ell)$$

die $k$-te symmetrische Potenz. Der Darstellungsraum $\mathrm{Sym}^k(V)$ (für $\dim V = 2$)
hat Dimension $k+1$, mit Basis $\{v_1^j \otimes v_2^{k-j} : j = 0, \ldots, k\}$.

### Frobenius-Spur

Seien $\alpha, \beta$ die Frobenius-Eigenwerte bei $p$. Dann:

$$\mathrm{tr}\bigl(\mathrm{Sym}^k(\rho)(\mathrm{Frob}_p)\bigr) = \sum_{j=0}^{k} \alpha^j \beta^{k-j}$$

Falls $\alpha \neq \beta$: $= \dfrac{\alpha^{k+1} - \beta^{k+1}}{\alpha - \beta}$

### Langlands-Funktorialität (Sym^k-Lifting)

Die Vermutung, dass $\mathrm{Sym}^k(\pi_E)$ eine automorphe Form auf $\mathrm{GL}_{k+1}$ ist,
ist für $k \leq 4$ bewiesen (Kim-Shahidi, 2002) und allgemein offen (aktives Forschungsgebiet).

---

## 6. Langlands-Korrespondenz `LanglandsCorrespondence`

### Das Langlands-Programm

Robert Langlands formulierte 1967 eine weitreichende Reihe von Vermutungen, die eine
"große vereinheitlichende Theorie" der Mathematik darstellen:

**Lokale Langlands-Korrespondenz** (bewiesen für $\mathrm{GL}_n$, Harris-Taylor 2001, Henniart 2000):
Isomorphismus zwischen $n$-dimensionalen Weil-Deligne-Darstellungen von $W_{\mathbb{Q}_p}$ und
irreduziblen glatten Darstellungen von $\mathrm{GL}_n(\mathbb{Q}_p)$.

**Globale Langlands-Korrespondenz** (weitgehend offen):
Für jede $n$-dimensionale motivische Galois-Darstellung $\rho: G_\mathbb{Q} \to \mathrm{GL}_n(\mathbb{C})$
existiert eine automorphe Darstellung $\pi$ von $\mathrm{GL}_n(\mathbb{A}_\mathbb{Q})$ mit

$$L(\rho, s) = L(\pi, s)$$

### Modularitätssatz (Wiles 1995)

Für jede semistabile elliptische Kurve $E/\mathbb{Q}$ existiert eine Modulform $f_E$ der Gewicht 2
mit $a_p(f_E) = a_p(E)$ für fast alle $p$. (Allgemein: Breuil-Conrad-Diamond-Taylor 2001.)

Dies ist der erste große Satz des Langlands-Programms für $\mathrm{GL}_2$.

### Ramanujan-Vermutung

Für den Tate-Modul einer elliptischen Kurve:

$$|a_p| \leq 2\sqrt{p} \quad \text{für gute Primzahlen } p$$

Allgemein für automorphe Formen: $|a_p| \leq 2 p^{(k-1)/2}$ (offen für allgemeine $k$).

### Funktionalgleichung

$$\Lambda(E, s) = \varepsilon_E \cdot \Lambda(E, 2-s)$$

wobei $\Lambda(E, s) = \left(\frac{\sqrt{N}}{2\pi}\right)^s \Gamma(s) L(E, s)$ und $\varepsilon_E \in \{+1, -1\}$.

---

## 7. Hilfsfunktion `galois_group_order`

| Polynom | Galois-Gruppe | Ordnung |
|---------|--------------|---------|
| $x^2 - 2$ | $\mathbb{Z}/2\mathbb{Z}$ | 2 |
| $x^3 - 2$ | $S_3$ | 6 |
| $x^4 - 2$ | $D_4$ (Diedergruppe) | 8 |
| $x^5 - x - 1$ | $S_5$ | 120 |
| $x^n - p$ ($p$ prim) | $\mathbb{Z}/n\mathbb{Z} \rtimes (\mathbb{Z}/n\mathbb{Z})^*$ | $n \cdot \varphi(n)$ |

Die Galois-Gruppe $\mathrm{Gal}(K/\mathbb{Q})$ wirkt treu und transitiv auf den Wurzeln
von $f$ und ist damit eine transitive Untergruppe der $S_n$ ($n = \deg f$).

---

## 8. Artin-Leiter `artin_conductor`

Der **Artin-Leiter** $N(\rho)$ misst die Verzweigung von $\rho$:

$$N(\rho) = \prod_p p^{f(\rho, p)}$$

mit lokalem Exponent

$$f(\rho, p) = \dim_{\mathbb{C}}(\rho/\rho^{I_p}) + \mathrm{Swan}(\rho, p)$$

- $I_p$ = Trägheitsgruppe (zahme Verzweigung: $f(\rho,p) \geq 1$; Verzweigung verschwindet: $f = 0$)
- $\mathrm{Swan}(\rho, p)$ = wilder Verzweigungsterm (nur für $p = \ell$ relevant)

---

## 9. Weil-Deligne-Darstellung `weil_deligne_representation_demo`

### Definition

Eine **Weil-Deligne-Darstellung** von $\mathbb{Q}_p$ ist ein Paar $(r, N)$:
- $r: W_p \to \mathrm{GL}_n(\mathbb{C})$ — Darstellung der Weil-Gruppe $W_p \subset \mathrm{Gal}(\bar{\mathbb{Q}}_p/\mathbb{Q}_p)$
- $N$ — nilpotente Matrix mit $r(\sigma) N r(\sigma)^{-1} = \|\sigma\|_p \cdot N$

### Reduktionstypen elliptischer Kurven

1. **Gute Reduktion** bei $p$: $N = 0$, Frobenius hat Eigenwerte $\alpha, \beta$ mit $|\alpha| = |\beta| = \sqrt{p}$
2. **Multiplikative Reduktion** (Tate-Kurve): $N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} \neq 0$
3. **Additive Reduktion**: Wilder Verzweigungsterm, komplizierter

---

## Literatur

- J.P. Serre: *Abelian ℓ-adic representations and elliptic curves* (1968)
- J.P. Serre: *Local Fields* (1979)
- J. Silverman: *The Arithmetic of Elliptic Curves* (1986)
- J. Silverman: *Advanced Topics in the Arithmetic of Elliptic Curves* (1994)
- R.P. Langlands: *Problems in the theory of automorphic forms* (1967)
- A. Wiles: *Modular elliptic curves and Fermat's last theorem* (1995)
- M. Harris, R. Taylor: *The geometry and cohomology of some simple Shimura varieties* (2001)
- P. Deligne: *La conjecture de Weil II* (1980)
