# cohen_lenstra.py – Dokumentation

**Modul:** `src/cohen_lenstra.py`
**Autor:** Michael Fuhrmann
**Stand:** Build 122 (2026-03-12)

---

## Übersicht

Das Modul implementiert die **Cohen-Lenstra-Heuristiken** (1984) für die statistische
Verteilung der **Idealklassengruppen** imaginär-quadratischer Zahlkörper
$K = \mathbb{Q}(\sqrt{-d})$.

---

## Mathematischer Hintergrund

### Klassengruppe und Klassenzahl

Für einen imaginär-quadratischen Zahlkörper $K = \mathbb{Q}(\sqrt{-d})$ ($d > 0$ quadratfrei)
ist die **Idealklassengruppe** $\text{Cl}(O_K)$ eine endliche abelsche Gruppe.

Ihre Ordnung $h(-d) = |\text{Cl}(O_K)|$ heißt **Klassenzahl**. Es gilt:
$$h(-d) = 1 \iff O_K \text{ ist ein Hauptidealbereich (PID)}$$

### Diskriminante

Die Diskriminante von $\mathbb{Q}(\sqrt{-d})$ ist:
$$D = \begin{cases} -d & \text{falls } d \equiv 3 \pmod{4} \\ -4d & \text{sonst} \end{cases}$$

### Dirichlet-Klassenanzahlformel

$$h(-d) = \frac{w \sqrt{|D|}}{2\pi} \cdot L(1, \chi_D)$$

wobei:
- $w = 2$ (für $d > 4$), $w = 4$ (für $d=1$), $w = 6$ (für $d=3$)
- $L(1, \chi_D) = \sum_{n=1}^{\infty} \chi_D(n)/n$ mit dem Kronecker-Symbol

### Heegner-Zahlen

Die einzigen squarefreen $d$ mit $h(-d) = 1$ sind die **9 Heegner-Zahlen**:
$$d \in \{1, 2, 3, 7, 11, 19, 43, 67, 163\}$$
*(Baker-Heegner-Stark-Theorem, 1966/1967)*

### Cohen-Lenstra-Heuristiken

Cohen und Lenstra (1984) postulierten:

**Gewichtung:** $\Pr[\text{Cl}(K) \cong A] \propto \frac{1}{|\text{Aut}(A)|}$

**p-Anteile:** Für ungerade Primzahlen $p$ und zufälliges imaginär-quadratisches $D$:

$$\Pr[p \mid h(D)] = 1 - \prod_{k=1}^{\infty} \left(1 - \frac{1}{p^k}\right)$$

**Numerische Werte:**

| $p$ | $\Pr[p \mid h]$ (theoretisch) |
|-----|-------------------------------|
| $2$ | $\approx 0.7112$ |
| $3$ | $\approx 0.4399$ |
| $5$ | $\approx 0.2386$ |
| $7$ | $\approx 0.1433$ |

*Hinweis:* Der 2-Anteil wird stärker durch die Genus-Theorie beeinflusst.

### Genus-Theorie

Die **Genus-Theorie** (Gauss) liefert eine untere Schranke:
$$h(-d) \geq 2^{t-1}$$
wobei $t$ = Anzahl der verschiedenen Primfaktoren von $D$.

Insbesondere: $h(-d)$ ist gerade, wenn $d$ mehr als einen Primfaktor hat.

### Bhargava-Shankar (2015)

Bhargava und Shankar bewiesen:
- Durchschnittlicher Rang elliptischer Kurven über $\mathbb{Q}$: $\text{avg}(\text{rank}) < 0.885$
- Verknüpfung mit Cohen-Lenstra: Selmer-Gruppen haben ähnliche Verteilung wie Klassengruppen

---

## Klassen

### `ClassNumberComputation`

Berechnet $h(-d)$ für imaginär-quadratische Körper.

| Methode | Beschreibung |
|---------|-------------|
| `class_number(d)` | Berechnet $h(-d)$ exakt (via reduzierte Formen) |
| `class_number_dirichlet_formula(d)` | Approximiert $h(-d)$ via Dirichlet-Formel |
| `batch_class_numbers(d_max)` | Berechnet alle $h(-d)$ für $d \leq d_{\max}$ |
| `genus_theory_lower_bound(d)` | Untere Schranke $2^{t-1}$ für $h(-d)$ |
| `heegner_numbers()` | Gibt die 9 Heegner-Zahlen zurück |

**Algorithmus für `class_number(d)`:**

Zählt reduzierte primitive binäre quadratische Formen $(a,b,c)$ mit Diskriminante $D$:
- Reduziertheitsbedingungen: $|b| \leq a \leq c$, $b \geq 0$ falls $|b|=a$ oder $a=c$
- Primitivität: $\gcd(a,b,c) = 1$

### `CohenLenstraHeuristics`

Implementiert die statistischen Vorhersagen von Cohen-Lenstra.

| Methode | Beschreibung |
|---------|-------------|
| `cohen_lenstra_probability(p, num_terms)` | Theoretische $\Pr[p \mid h(-d)]$ |
| `empirical_probability(p)` | Empirische Häufigkeit aus Daten |
| `compare_predictions(primes)` | Vergleich Theorie vs. Empirie |
| `cohen_lenstra_weight(group_order)` | Gewicht $1/\varphi(n)$ für $\mathbb{Z}/n\mathbb{Z}$ |
| `verify_heegner_numbers()` | Verifiziert alle 9 Heegner-Zahlen |
| `class_number_distribution()` | Statistische Auswertung der $h(-d)$ |
| `divisibility_statistics(primes)` | Teilbarkeitshäufigkeiten |
| `bhargava_shankar_summary()` | Zusammenfassung Bhargava-Shankar |
| `p_rank_distribution(p)` | Theoretische $p$-Rang-Verteilung |

---

## Verwendungsbeispiel

```python
from cohen_lenstra import ClassNumberComputation, CohenLenstraHeuristics

# Klassenzahl berechnen
cnc = ClassNumberComputation()
print(cnc.class_number(5))    # h(-5) = 2
print(cnc.class_number(163))  # h(-163) = 1 (Heegner-Zahl)
print(cnc.heegner_numbers())  # [1, 2, 3, 7, 11, 19, 43, 67, 163]

# Cohen-Lenstra-Vorhersage
cl = CohenLenstraHeuristics(d_max=1000)
prob3 = cl.cohen_lenstra_probability(3)
print(f'Pr[3|h(-d)] ≈ {prob3:.4f}')  # ≈ 0.4399

# Vergleich Theorie vs. Empirie
for entry in cl.compare_predictions([3, 5]):
    print(f"p={entry['p']}: theoret.={entry['theoretical']:.4f}, "
          f"empir.={entry['empirical']:.4f}")
```

---

## Bekannte Einschränkungen

- `class_number(d)` verwendet squarefree kernel von $d$; für nicht-squarefrees $d$ wird der kernel benutzt
- Empirische Konvergenz von Cohen-Lenstra erfordert große $d_{\max}$ (> 1000)
- Die Dirichlet-Formel-Approximation nutzt 10000 Summanden (hohe Genauigkeit für kleine $d$)

---

## Referenzen

1. Cohen, H. & Lenstra, H.W. (1984). *Heuristics on class groups of number fields.* LNM **1068**, 33–62.
2. Baker, A. (1966). *Linear forms in the logarithms of algebraic numbers.* Mathematika **13**, 204–216.
3. Stark, H.M. (1967). *A complete determination of the complex quadratic fields of class-number one.* Mich. Math. J. **14**, 1–27.
4. Bhargava, M. & Shankar, A. (2015). *Ternary cubic forms having bounded invariants.* Ann. Math. **181**, 587–621.
5. Cohen, H. (1993). *A Course in Computational Algebraic Number Theory.* Springer.
