# Modul: Operations Research (`operations_research.py`)

> **Autor:** Kurt Ingwer | **Letzte Änderung:** 2026-03-10

## Überblick

Das Modul `operations_research.py` implementiert zentrale Verfahren der Operationsforschung in sieben Klassen. Es deckt lineare und ganzzahlige Optimierung (Simplex, Branch-and-Bound, Rucksackproblem, TSP), Warteschlangentheorie (M/M/1, M/M/s, M/M/1/K, Erlang), Markov-Entscheidungsprozesse (Wertiteration, Politikiteration), Netzwerkflüsse (Ford-Fulkerson, Dijkstra, Kruskal), Spieltheorie (Nash-Gleichgewicht, Minimax) sowie stochastische Optimierung (SGD, Simulated Annealing) ab.

---

## Mathematischer Hintergrund

### Lineare Optimierung (LP)

**Standardform** (Minimierung):
$$\min_{\mathbf{x}} \; \mathbf{c}^T \mathbf{x} \quad \text{s.t.} \quad A\mathbf{x} \leq \mathbf{b},\; \mathbf{x} \geq \mathbf{0}$$

**Simplex-Methode:** Traversiert die Ecken des konvexen Polyeders (zulässige Menge) in Richtung fallender Zielfunktion. Da LP-Optima in Ecken liegen, ist Simplex exakt (aber im Worst-Case exponentiell; in der Praxis sehr effizient).

**Duales Problem:**
$$\max_{\mathbf{y}} \; \mathbf{b}^T \mathbf{y} \quad \text{s.t.} \quad A^T \mathbf{y} \leq \mathbf{c},\; \mathbf{y} \geq \mathbf{0}$$

**Starker Dualitätssatz:** Primal und Dual haben den selben Optimalwert (falls beide zulässig).

**Sensitivitätsanalyse:** Wie ändert sich die optimale Lösung bei kleinen Perturbationen von $\mathbf{c}$ oder $\mathbf{b}$? (Schattenpreise = Dualvariablen)

**Transportproblem:** Spezielle LP-Struktur – $m$ Angebotsorte, $n$ Nachfrageorte, Minimierung der Transportkosten.

**Zuweisungsproblem:** $n \times n$-Kostematrix, bijektive Zuweisung minimaler Kosten. Optimal per ungarischem Algorithmus ($O(n^3)$).

### Ganzzahlige Optimierung (IP)

**Branch-and-Bound:** Löse LP-Relaxierung; verzweige an Variablen mit gebrochenem Wert.

**Gomory-Schnittebenen:** Füge gültige Ungleichungen hinzu, die alle ganzzahligen Punkte enthalten, aber den LP-Optimalpunkt ausschließen.

**0-1-Rucksack:**
$$\max \sum_{i=1}^n v_i x_i \quad \text{s.t.} \quad \sum_{i=1}^n w_i x_i \leq W,\; x_i \in \{0,1\}$$

Optimal per DP in $O(nW)$; NP-vollständig im allgemeinen Sinne.

**Fraktionaler Rucksack:** $x_i \in [0,1]$ – lösbar durch Greedy nach Wert/Gewicht-Verhältnis.

**TSP** (Handlungsreisendenproblem): NP-schwer. Nearest-Neighbor-Heuristik: $O(n^2)$, Approximationsfaktor unbeschränkt. 2-Opt-Verbesserung: vertausche Kanten, solange die Tour kürzer wird.

### Warteschlangentheorie

Notation: $\lambda$ = Ankunftsrate, $\mu$ = Servicerate, $\rho = \lambda/\mu$ (Auslastung), $s$ = Anzahl Server.

**M/M/1-Warteschlange** ($s=1$, $\rho < 1$):
$$L = \frac{\rho}{1-\rho}, \quad W = \frac{1}{\mu - \lambda}, \quad P_0 = 1 - \rho$$

**M/M/s-Warteschlange** ($s$ Server):
$$P_0 = \left[\sum_{n=0}^{s-1}\frac{(\lambda/\mu)^n}{n!} + \frac{(\lambda/\mu)^s}{s!(1-\rho/s)}\right]^{-1}$$

**M/M/1/K** (endliche Kapazität $K$):
$$P_n = \frac{(1-\rho)\rho^n}{1-\rho^{K+1}}, \quad n = 0,1,\ldots,K$$

**Erlang-C-Formel** (Wahrscheinlichkeit Warten in M/M/s):
$$C(s,\lambda/\mu) = \frac{\frac{(\lambda/\mu)^s}{s!} \cdot \frac{s}{s-\lambda/\mu}}{\sum_{n=0}^{s-1}\frac{(\lambda/\mu)^n}{n!} + \frac{(\lambda/\mu)^s}{s!} \cdot \frac{s}{s-\lambda/\mu}}$$

**Little's Law:** $L = \lambda W$ (gilt für stationäre Systeme sehr allgemein).

### Markov-Entscheidungsprozesse (MDP)

**Bellman-Gleichung** (optimaler Wert):
$$V^*(s) = \max_{a \in A} \left[ R(s,a) + \gamma \sum_{s'} P(s'|s,a) V^*(s') \right]$$

mit Diskontierungsfaktor $\gamma \in [0,1)$.

**Wertiteration:** Iteriere $V_{k+1}(s) = \max_a[R(s,a) + \gamma \sum_{s'} P V_k(s')]$ bis Konvergenz ($O(|S|^2|A|/\varepsilon)$).

**Politikiteration:** Werte aktuelle Politik aus (löse lineares System), dann verbessere die Politik. Konvergiert in endlich vielen Schritten.

### Netzwerkflüsse

**Max-Flow** (Ford-Fulkerson): Finde augmentierende Pfade per BFS (Edmonds-Karp, $O(VE^2)$) und erhöhe den Fluss.

**Min-Cut** (Max-Flow-Min-Cut-Theorem): Maximaler Fluss = minimale Schnittkapazität.

**Kürzeste Wege:**
- **Dijkstra:** $O((V+E)\log V)$ mit Prioritätswarteschlange; nur positive Kantengewichte.
- **Bellman-Ford:** $O(VE)$; erkennt negative Zyklen.

**Minimaler Spannbaum:**
- **Kruskal:** Sortiere Kanten nach Gewicht, füge greedy hinzu (Union-Find), $O(E\log E)$.
- **Prim:** Wachse vom Startknoten aus, $O((V+E)\log V)$.

### Spieltheorie

**Nash-Gleichgewicht** für 2×2-Spiele: Strategieprofil $(x^*, y^*)$, bei dem kein Spieler durch einseitige Abweichung gewinnt. Gemischte Strategie: $x^* = (d-b)/[(a-c)+(d-b)]$.

**Minimax-Theorem** (von Neumann 1928): Für Zwei-Personen-Nullsummenspiele mit Auszahlungsmatrix $A$:
$$\max_x \min_y x^T A y = \min_y \max_x x^T A y$$

**Gefangenendilemma:** Klassisches Beispiel für Konflikt zwischen individuellem und kollektivem Optimum.

**Pareto-Optimalität:** Ein Ergebnis ist Pareto-optimal, wenn kein Spieler besser gestellt werden kann, ohne einen anderen schlechter zu stellen.

### Stochastische Optimierung

**Stochastischer Gradientenabstieg (SGD):**
$$\theta_{k+1} = \theta_k - \eta_k \nabla_\theta f(\theta_k, \xi_k)$$

**Simulated Annealing:** Metropolis-Heuristik – akzeptiere schlechtere Lösungen mit Wahrscheinlichkeit $e^{-\Delta E/T}$, kühle $T$ langsam ab.

**Zeitungshändlerproblem (Newsvendor):**
$$Q^* = F^{-1}\!\left(\frac{c_u}{c_u + c_o}\right)$$

mit $c_u$ = Unterkosten (entgangener Gewinn), $c_o$ = Überkosten (Lagerkosten).

---

## Klassen und Methoden

### `LinearProgramming`

Lineare Optimierung, Transportprobleme, Zuweisungsprobleme.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `simplex_method` | `(c, A_ub, b_ub, ...) → dict` | LP-Minimierung via `scipy.optimize.linprog` |
| `dual_problem` | `(c, A_ub, b_ub) → dict` | Duales LP (Maximierung) |
| `sensitivity_analysis` | `(c, A_ub, b_ub, ...) → dict` | Sensitivitätsanalyse (Schattenpreise) |
| `transportation_problem` | `(supply, demand, costs) → dict` | Transportproblem (LP-Formulierung) |
| `assignment_problem` | `(cost_matrix) → dict` | Zuweisungsproblem (Ungarischer Algorithmus) |
| `two_phase_simplex_demo` | `(c, A_eq, b_eq, ...) → dict` | Zwei-Phasen-Simplex (Demo) |

### `IntegerProgramming`

Ganzzahlige Optimierung und kombinatorische Algorithmen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `branch_and_bound_demo` | `(c, A_ub, b_ub, ...) → dict` | Branch-and-Bound (Demo) |
| `cutting_plane_demo` | `() → dict` | Gomory-Schnittebenen (Demo) |
| `knapsack_01` | `(values, weights, capacity) → dict` | 0-1-Rucksack per DP |
| `knapsack_fractional` | `(values, weights, capacity) → dict` | Fraktionaler Rucksack (Greedy) |
| `traveling_salesman_nearest_neighbor` | `(dist_matrix, start=0) → dict` | TSP Nearest-Neighbor-Heuristik |
| `traveling_salesman_2opt` | `(dist_matrix, initial_tour) → dict` | TSP 2-Opt-Verbesserung |

### `QueueingTheory`

Stochastische Warteschlangenmodelle.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `mm1_queue` | `(arrival_rate, service_rate) → dict` | M/M/1: $L, W, L_q, W_q, P_0$ |
| `mms_queue` | `(arrival_rate, service_rate, s) → dict` | M/M/s: Leistungsmaße mit $s$ Servern |
| `mm1k_queue` | `(arrival_rate, service_rate, K) → dict` | M/M/1/K: endliche Kapazität |
| `erlang_c` | `(s, lam_over_mu) → dict` | Erlang-C-Formel (Warte-Wahrscheinlichkeit) |
| `erlang_b` | `(s, lam_over_mu) → float` | Erlang-B-Formel (Blockierungs-Wahrscheinlichkeit) |
| `little_law_check` | `(L, lambda_, W) → dict` | Verifikation von $L = \lambda W$ |

### `MarkovDecisionProcess`

Markov-Entscheidungsprozesse und Verstärkungslernen-Grundlagen.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `bellman_equation` | `(V, s) → float` | Bellman-Optimalitätsgleichung für Zustand $s$ |
| `value_iteration` | `(max_iter, tol) → dict` | Wertiteration bis Konvergenz |
| `policy_iteration` | `(max_iter) → dict` | Politikiteration (exakte Politikauswertung) |
| `optimal_policy` | `() → dict` | Optimale Politik aus konvergiertem $V^*$ |

### `NetworkFlows`

Netzwerkfluss-Algorithmen und Graphoptimierung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `max_flow_ford_fulkerson` | `(graph, source, sink) → dict` | Max-Flow (Edmonds-Karp via BFS) |
| `min_cut` | `(graph, source, sink) → dict` | Min-Cut aus Max-Flow-Ergebnis |
| `shortest_path_dijkstra` | `(graph, source) → dict` | Dijkstra ($O((V+E)\log V)$) |
| `shortest_path_bellman_ford` | `(graph, source) → dict` | Bellman-Ford, erkennt negative Zyklen |
| `minimum_spanning_tree_kruskal` | `(n, edges) → dict` | Kruskal (Union-Find, $O(E\log E)$) |
| `minimum_spanning_tree_prim` | `(adj_matrix) → dict` | Prim (Prioritätswarteschlange) |

### `GameTheory`

Spieltheorie – Nash-Gleichgewichte, Minimax, Strategien.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `nash_equilibrium_2x2` | `(payoff_A, payoff_B) → dict` | Nash-GGW für 2×2-Bimatrixspiele |
| `minimax_theorem` | `(A) → dict` | Minimax-Wert für Nullsummenspiel |
| `prisoners_dilemma` | `() → dict` | Gefangenendilemma – Analyse |
| `dominant_strategy` | `(payoff_matrix) → dict` | Dominante Strategie (iterierte Eliminierung) |
| `pareto_optimal` | `(outcomes) → dict` | Pareto-optimale Ergebnisse |

### `StochasticOptimization`

Stochastische und metaheuristische Optimierung.

| Methode | Signatur | Beschreibung |
|---------|----------|--------------|
| `stochastic_gradient_descent` | `(f, grad_f, x0, n_iter, lr, ...) → dict` | SGD mit Lernratenplan |
| `simulated_annealing` | `(f, x0, T_init, T_min, alpha, n_iter) → dict` | Simulated Annealing |
| `newsvendor_problem` | `(demand_dist, c_underage, c_overage, ...) → dict` | Zeitungshändlerproblem: $Q^*$ |
| `stochastic_programming_demo` | `() → dict` | Zwei-Stufen-Stochastische Programmierung |

---

## Beispiele

```python
from operations_research import (LinearProgramming, IntegerProgramming,
                                  QueueingTheory, MarkovDecisionProcess,
                                  NetworkFlows, GameTheory)

# -- LP: Produktionsplanung --
lp = LinearProgramming()
c = [-3, -5]                          # Maximiere 3x+5y → Minimiere -3x-5y
A = [[1, 0], [0, 2], [3, 2]]
b = [4, 12, 18]
result = lp.simplex_method(c, A, b)
print(result["x"], result["obj"])     # [2, 6], obj = -36

# -- 0-1-Rucksack --
ip = IntegerProgramming()
res = ip.knapsack_01(values=[60, 100, 120], weights=[10, 20, 30], capacity=50)
print(res["selected_items"], res["total_value"])  # [1,2], 220

# -- M/M/1-Warteschlange --
qt = QueueingTheory()
stats = qt.mm1_queue(arrival_rate=3.0, service_rate=5.0)
print(stats["L"], stats["W"])         # L = 1.5, W = 0.5

# -- MDP: Wertiteration --
mdp = MarkovDecisionProcess(
    states=[0,1,2], actions=["left","right"],
    transition={(0,"right"):[(1,1.0)], (1,"right"):[(2,1.0)], (2,"right"):[(2,1.0)],
                (0,"left"):[(0,1.0)], (1,"left"):[(0,1.0)], (2,"left"):[(1,1.0)]},
    reward={(0,"right"):0, (1,"right"):0, (2,"right"):10,
            (0,"left"):0,  (1,"left"):0,  (2,"left"):0},
    gamma=0.9
)
vi = mdp.value_iteration()
print(vi["V"])

# -- Dijkstra --
nf = NetworkFlows()
graph = {0:{1:4,2:1}, 1:{3:1}, 2:{1:2,3:5}, 3:{}}
res = nf.shortest_path_dijkstra(graph, source=0)
print(res["distances"])               # {0:0, 1:3, 2:1, 3:4}

# -- Nash-Gleichgewicht Gefangenendilemma --
gt = GameTheory()
pd = gt.prisoners_dilemma()
print(pd["nash_equilibrium"])         # (Defect, Defect)

# -- Simulated Annealing: Minimierung --
so = StochasticOptimization()
res = so.simulated_annealing(f=lambda x: (x-3)**2, x0=0.0, T_init=100, T_min=0.01, alpha=0.99)
print(f"x* ≈ {res['x_best']:.3f}")   # ≈ 3.0
```

---

## Tests

**Testdatei:** `tests/test_operations_research.py`
**Abdeckung:** LP-Optimalität, Dualitätslücke, Rucksack-DP-Korrektheit, TSP-Tourvalidierung, Warteschlangenformeln, MDP-Konvergenz, Max-Flow-Min-Cut, Dijkstra-Korrektheit, Nash-Gleichgewichte, Simulated-Annealing-Konvergenz

---

## Implementierungshinweise

- **Simplex:** Nutzt `scipy.optimize.linprog` (intern: Revised-Simplex oder HiGHS); eigene Simplex-Demo für Zwei-Phasen-Methode.
- **Ungarischer Algorithmus:** Nutzt `scipy.optimize.linear_sum_assignment` ($O(n^3)$).
- **Ford-Fulkerson:** Implementiert mit BFS-Pfadsuche (Edmonds-Karp-Variante) für polynomielle Laufzeit $O(VE^2)$.
- **Union-Find:** Mit Pfadkomprimierung und Union-by-Rank für effizientes Kruskal.
- **Bellman-Ford:** Erkennt negative Zyklen durch $|V|$-te Iteration.
- **Simulated Annealing:** Geometrischer Kühlplan $T_{k+1} = \alpha \cdot T_k$; bei diskreten Problemen (z.B. TSP) muss die Nachbarschaftsstruktur angepasst werden.
- **Newsvendor:** Kritisches Verhältnis $c_u/(c_u + c_o)$ entspricht dem Quantil der optimalen Bestellmenge bei bekannter Nachfrageverteilung.
