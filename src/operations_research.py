"""
Operations Research – Operationsforschung.

Implementiert zentrale Verfahren der Operationsforschung:
    - Lineare Optimierung (LP): Simplex, duales Problem, Sensitivitätsanalyse,
      Transportproblem, Zuweisungsproblem, Zwei-Phasen-Simplex (Demo)
    - Ganzzahlige Optimierung (IP): Branch-and-Bound (Demo), Gomory-Schnittebenen (Demo),
      0-1-Rucksackproblem, fraktionales Rucksackproblem, TSP Nearest-Neighbor, 2-Opt
    - Warteschlangentheorie: M/M/1, M/M/s, M/M/1/K, Erlang-C, Erlang-B, Little's Law
    - Markov-Entscheidungsprozesse (MDP): Wertiteration, Politikiteration, Bellman-Gleichung
    - Netzwerkflüsse: Max-Flow (Ford-Fulkerson), Min-Cut, Dijkstra, Bellman-Ford,
      Kruskal MST, Prim MST
    - Spieltheorie: Nash-Gleichgewicht 2×2, Minimax-Theorem, Gefangenendilemma,
      dominante Strategie, Pareto-Optimalität
    - Stochastische Optimierung: SGD, Simulated Annealing, Zeitungshändlerproblem,
      stochastische Programmierung

Mathematischer Hintergrund:
    Operations Research kombiniert mathematische Modellierung, Statistik und
    algorithmische Methoden, um Entscheidungsprozesse zu optimieren.

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-10
@lastModified: 2026-03-10
"""

import math
import heapq
import numpy as np
from collections import defaultdict, deque
from typing import Callable, Dict, List, Optional, Tuple, Union
from scipy.optimize import linprog
from scipy.optimize import linear_sum_assignment


# ===========================================================================
# 1. LINEARE OPTIMIERUNG
# ===========================================================================

class LinearProgramming:
    """
    Lineare Optimierung (LP).

    Löst Probleme der Form:
        min  cᵀx
        s.t. A_ub · x ≤ b_ub
             A_eq · x  = b_eq
             x ≥ 0

    Verwendet scipy.optimize.linprog als numerisches Backend.

    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """

    @staticmethod
    def simplex_method(
        c: List[float],
        A_ub: List[List[float]],
        b_ub: List[float],
        A_eq: Optional[List[List[float]]] = None,
        b_eq: Optional[List[float]] = None,
        bounds: Optional[List[Tuple]] = None
    ) -> Dict:
        """
        Simplex-Methode via scipy.optimize.linprog.

        Löst das LP-Problem:
            min  cᵀx
            s.t. A_ub · x ≤ b_ub
                 A_eq · x  = b_eq  (optional)
                 lb ≤ x ≤ ub       (Standard: x ≥ 0)

        @param c:      Koeffizientenvektor der Zielfunktion (Minimierung)
        @param A_ub:   Koeffizientenmatrix der Ungleichungsnebenbedingungen
        @param b_ub:   Rechte Seite der Ungleichungsnebenbedingungen
        @param A_eq:   Koeffizientenmatrix der Gleichungsnebenbedingungen (optional)
        @param b_eq:   Rechte Seite der Gleichungsnebenbedingungen (optional)
        @param bounds: Liste von (lb, ub) für jede Variable; None = (0, None)
        @return: Dict mit Schlüsseln: 'x', 'fun', 'status', 'message', 'success'
        @lastModified: 2026-03-10
        """
        # Standardschranken: alle Variablen ≥ 0
        if bounds is None:
            bounds = [(0, None)] * len(c)

        # Leere Listen als None behandeln (linprog erwartet None wenn keine Nebenbedingung)
        A_ub_arg = A_ub if (A_ub is not None and len(A_ub) > 0) else None
        b_ub_arg = b_ub if (b_ub is not None and len(b_ub) > 0) else None

        # scipy.linprog aufrufen (minimiert cᵀx)
        result = linprog(
            c=c,
            A_ub=A_ub_arg,
            b_ub=b_ub_arg,
            A_eq=A_eq,
            b_eq=b_eq,
            bounds=bounds,
            method='highs'
        )

        return {
            'x':       result.x,
            'fun':     result.fun,
            'status':  result.status,
            'message': result.message,
            'success': result.success
        }

    @staticmethod
    def dual_problem(
        c: List[float],
        A: List[List[float]],
        b: List[float]
    ) -> Dict:
        """
        Duales LP aufstellen.

        Primal (Minimierung):
            min  cᵀx
            s.t. A · x ≥ b
                 x ≥ 0

        Das duale Problem (Maximierung) lautet:
            max  bᵀy
            s.t. Aᵀ · y ≤ c
                 y ≥ 0

        Schwacher Dualitätssatz: bᵀy ≤ cᵀx für zulässige y, x.
        Starker Dualitätssatz: Optimalwert von Primal = Optimalwert von Dual.

        @param c: Zielfunktionskoeffizienten des Primal-Problems (Minimierung)
        @param A: Constraint-Matrix des Primal-Problems (A·x ≥ b)
        @param b: Rechte Seite des Primal-Problems
        @return: Dict mit Dual-Formulierung und Lösung
        @lastModified: 2026-03-10
        """
        A_np = np.array(A, dtype=float)
        b_np = np.array(b, dtype=float)
        c_np = np.array(c, dtype=float)

        m, n = A_np.shape  # m Nebenbedingungen, n Variablen

        # Dual: max bᵀy  ⟺  min -bᵀy
        # s.t. Aᵀy ≤ c  (y ≥ 0)
        c_dual   = (-b_np).tolist()          # Zielfunktion negiert (linprog minimiert)
        A_ub_dual = A_np.T.tolist()          # Aᵀ · y ≤ c
        b_ub_dual = c_np.tolist()
        bounds_dual = [(0, None)] * m

        result = linprog(
            c=c_dual,
            A_ub=A_ub_dual,
            b_ub=b_ub_dual,
            bounds=bounds_dual,
            method='highs'
        )

        return {
            'dual_variables': result.x,
            'dual_objective': -result.fun if result.fun is not None else None,
            'success':        result.success,
            'message':        result.message,
            'dual_c':         c_dual,
            'dual_A_ub':      A_ub_dual,
            'dual_b_ub':      b_ub_dual
        }

    @staticmethod
    def sensitivity_analysis(
        c: List[float],
        A_ub: List[List[float]],
        b_ub: List[float]
    ) -> Dict:
        """
        Sensitivitätsanalyse (Schattenpreise / Shadow Prices).

        Schattenpreise messen, wie sich der Optimalwert ändert, wenn die
        rechte Seite b einer Nebenbedingung um 1 erhöht wird.

        Numerisch approximiert via:
            shadow_price[i] = (f(b + δ·eᵢ) - f(b)) / δ

        @param c:     Zielfunktionskoeffizienten
        @param A_ub:  Ungleichungsmatrix
        @param b_ub:  Rechte Seite
        @return: Dict mit 'shadow_prices', 'optimal_value', 'optimal_x'
        @lastModified: 2026-03-10
        """
        b_np = np.array(b_ub, dtype=float)
        m    = len(b_ub)
        n    = len(c)
        bounds = [(0, None)] * n

        # Basis-LP lösen
        base = linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')
        if not base.success:
            return {
                'shadow_prices':  None,
                'optimal_value':  None,
                'optimal_x':      None,
                'message':        base.message
            }

        f0 = base.fun
        delta = 1e-4   # Perturbationsgröße

        shadow_prices = []
        for i in range(m):
            # b um delta in Richtung i erhöhen
            b_pert = b_np.copy()
            b_pert[i] += delta
            res = linprog(c=c, A_ub=A_ub, b_ub=b_pert.tolist(), bounds=bounds, method='highs')
            if res.success:
                # Schattenpreis = (f_neu - f_alt) / delta
                shadow_prices.append((res.fun - f0) / delta)
            else:
                shadow_prices.append(None)

        return {
            'shadow_prices': shadow_prices,
            'optimal_value': f0,
            'optimal_x':     base.x
        }

    @staticmethod
    def transportation_problem(
        supply: List[float],
        demand: List[float],
        costs: List[List[float]]
    ) -> Dict:
        """
        Transportproblem (als LP formuliert und gelöst).

        Minimiert Transportkosten:
            min  Σᵢ Σⱼ cᵢⱼ · xᵢⱼ
            s.t. Σⱼ xᵢⱼ = supplyᵢ  (Angebotsrestriktion)
                 Σᵢ xᵢⱼ = demandⱼ  (Nachfragerestriktion)
                 xᵢⱼ ≥ 0

        Bedingung: Σ supply = Σ demand (balanciertes Problem).
        Bei Unbalance wird ein fiktiver Lieferant/Empfänger hinzugefügt.

        @param supply: Angebotsmengen der m Lieferanten
        @param demand: Nachfragemengen der n Empfänger
        @param costs:  Kostenmatrix (m × n)
        @return: Dict mit 'flow_matrix', 'total_cost', 'success'
        @lastModified: 2026-03-10
        """
        supply = list(supply)
        demand = list(demand)
        costs  = [list(row) for row in costs]

        m = len(supply)
        n = len(demand)

        total_supply = sum(supply)
        total_demand = sum(demand)

        # Unbalanciertes Problem ausgleichen
        if abs(total_supply - total_demand) > 1e-9:
            if total_supply > total_demand:
                # Fiktiven Empfänger mit 0-Kosten hinzufügen
                diff = total_supply - total_demand
                demand.append(diff)
                for row in costs:
                    row.append(0.0)
                n += 1
            else:
                # Fiktiven Lieferanten mit 0-Kosten hinzufügen
                diff = total_demand - total_supply
                supply.append(diff)
                costs.append([0.0] * n)
                m += 1

        # LP-Formulierung: Variablenvektor x = [x_00, x_01, ..., x_{m-1,n-1}]
        # Kosten flachgedrückt
        c_lp = [costs[i][j] for i in range(m) for j in range(n)]

        # Gleichungsnebenbedingungen: Angebot- und Nachfragerestriktionen
        A_eq = []
        b_eq = []

        # Angebotsrestriktionen: Σⱼ x_{ij} = supply_i
        for i in range(m):
            row = [0.0] * (m * n)
            for j in range(n):
                row[i * n + j] = 1.0
            A_eq.append(row)
            b_eq.append(supply[i])

        # Nachfragerestriktionen: Σᵢ x_{ij} = demand_j
        for j in range(n):
            row = [0.0] * (m * n)
            for i in range(m):
                row[i * n + j] = 1.0
            A_eq.append(row)
            b_eq.append(demand[j])

        bounds = [(0, None)] * (m * n)
        result = linprog(c=c_lp, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

        if result.success:
            flow_matrix = np.array(result.x).reshape(m, n)
            return {
                'flow_matrix': flow_matrix,
                'total_cost':  result.fun,
                'success':     True
            }
        else:
            return {
                'flow_matrix': None,
                'total_cost':  None,
                'success':     False,
                'message':     result.message
            }

    @staticmethod
    def assignment_problem(cost_matrix: List[List[float]]) -> Dict:
        """
        Zuweisungsproblem (Ungarische Methode via scipy).

        Findet eine bijektive Zuweisung von m Agenten zu n Aufgaben,
        die die Gesamtkosten minimiert:
            min  Σᵢ cᵢ,σ(ᵢ)  über alle Permutationen σ

        Nutzt den Kuhn-Munkres-Algorithmus (Ungarische Methode) via
        scipy.optimize.linear_sum_assignment (O(n³)).

        @param cost_matrix: Kostenmatrix (m × n), quadratisch oder rechteckig
        @return: Dict mit 'row_ind', 'col_ind', 'total_cost'
        @raises ValueError: Wenn cost_matrix leer ist
        @lastModified: 2026-03-10
        """
        C = np.array(cost_matrix, dtype=float)
        if C.size == 0:
            raise ValueError("cost_matrix darf nicht leer sein")

        # Ungarische Methode via scipy
        row_ind, col_ind = linear_sum_assignment(C)

        total_cost = C[row_ind, col_ind].sum()

        return {
            'row_ind':    row_ind.tolist(),
            'col_ind':    col_ind.tolist(),
            'total_cost': float(total_cost)
        }

    @staticmethod
    def two_phase_simplex_demo(
        c: List[float],
        A: List[List[float]],
        b: List[float]
    ) -> Dict:
        """
        Zwei-Phasen-Simplex-Methode (Demo).

        Phase 1: Finde eine zulässige Startlösung, indem künstliche Variablen
                 eingeführt und ihre Summe minimiert wird.
        Phase 2: Optimiere die ursprüngliche Zielfunktion, beginnend mit der
                 zulässigen Startlösung aus Phase 1.

        Formulierung: min cᵀx  s.t. A·x = b, x ≥ 0.

        @param c: Zielfunktionskoeffizienten
        @param A: Gleichungsmatrix (A·x = b)
        @param b: Rechte Seite (alle bᵢ ≥ 0 vorausgesetzt)
        @return: Dict mit 'phase1_result', 'phase2_result', 'feasible', 'optimal_value'
        @lastModified: 2026-03-10
        """
        m = len(b)   # Anzahl Zeilen
        n = len(c)   # Anzahl Originalvariablen

        # --- Phase 1 ---
        # Minimiere Σ künstliche Variablen aᵢ
        # s.t. A·x + I·a = b, x ≥ 0, a ≥ 0
        c_phase1 = [0.0] * n + [1.0] * m
        A_eq_phase1 = []
        for i, row in enumerate(A):
            # Originalvariablen + Einheitsmatrix für künstliche Variablen
            art_row = [0.0] * m
            art_row[i] = 1.0
            A_eq_phase1.append(list(row) + art_row)

        bounds_phase1 = [(0, None)] * (n + m)
        res1 = linprog(
            c=c_phase1, A_eq=A_eq_phase1, b_eq=b,
            bounds=bounds_phase1, method='highs'
        )

        feasible = res1.success and abs(res1.fun) < 1e-8

        if not feasible:
            return {
                'phase1_result': res1.message,
                'phase2_result': None,
                'feasible':      False,
                'optimal_value': None
            }

        # --- Phase 2 ---
        # Minimiere cᵀx  s.t. A·x = b, x ≥ 0
        bounds_phase2 = [(0, None)] * n
        res2 = linprog(
            c=c, A_eq=A, b_eq=b,
            bounds=bounds_phase2, method='highs'
        )

        return {
            'phase1_result': 'Zulässige Startlösung gefunden',
            'phase2_result': res2.message,
            'feasible':      True,
            'optimal_value': res2.fun if res2.success else None,
            'optimal_x':     res2.x if res2.success else None
        }


# ===========================================================================
# 2. GANZZAHLIGE OPTIMIERUNG
# ===========================================================================

class IntegerProgramming:
    """
    Ganzzahlige Optimierung (IP).

    Methoden für kombinatorische Optimierungsprobleme:
    Branch-and-Bound, Gomory-Schnittebenen, Rucksackproblem, TSP-Heuristiken.

    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """

    @staticmethod
    def branch_and_bound_demo(
        c: List[float],
        A_ub: List[List[float]],
        b_ub: List[float],
        bounds: Optional[List[Tuple]] = None
    ) -> Dict:
        """
        Branch-and-Bound für ganzzahlige LP (Demo, vereinfacht).

        Strategie: Best-First-Search mit LP-Relaxation als Schranke.
        Verzweigungsregel: Wähle erste nicht-ganzzahlige Variable.

        @param c:      Zielfunktionskoeffizienten (Minimierung)
        @param A_ub:   Ungleichungsmatrix
        @param b_ub:   Rechte Seite
        @param bounds: Variablenschranken (Standard: 0 ≤ x ≤ ∞)
        @return: Dict mit 'optimal_x', 'optimal_value', 'iterations', 'success'
        @lastModified: 2026-03-10
        """
        n = len(c)
        if bounds is None:
            bounds = [(0, None)] * n

        # LP-Relaxation als Ausgangspunkt
        best_value  = math.inf
        best_x      = None
        iterations  = 0

        # Stapel für Branch-and-Bound-Knoten: (bounds_für_diesen_Knoten,)
        stack = [list(bounds)]

        while stack:
            current_bounds = stack.pop()
            iterations += 1

            # LP-Relaxation lösen
            res = linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=current_bounds, method='highs')

            if not res.success:
                # Unzulässiger Teilbaum → überspringen
                continue

            if res.fun >= best_value:
                # Schlechter als aktuelle beste Lösung → Pruning
                continue

            # Prüfen ob alle Variablen ganzzahlig
            x = res.x
            branch_var = -1
            for i in range(n):
                if abs(x[i] - round(x[i])) > 1e-6:
                    branch_var = i
                    break

            if branch_var == -1:
                # Alle Variablen ganzzahlig → neue beste Lösung
                best_value = res.fun
                best_x     = [round(xi) for xi in x]
            else:
                # Verzweigung an branch_var
                frac_val = x[branch_var]

                # Linker Ast: x[branch_var] ≤ floor(frac_val)
                bounds_left = list(current_bounds)
                lb, ub = bounds_left[branch_var]
                bounds_left[branch_var] = (lb, math.floor(frac_val))
                stack.append(bounds_left)

                # Rechter Ast: x[branch_var] ≥ ceil(frac_val)
                bounds_right = list(current_bounds)
                bounds_right[branch_var] = (math.ceil(frac_val), ub)
                stack.append(bounds_right)

        return {
            'optimal_x':     best_x,
            'optimal_value': best_value if best_x is not None else None,
            'iterations':    iterations,
            'success':       best_x is not None
        }

    @staticmethod
    def cutting_plane_demo() -> Dict:
        """
        Gomory-Schnittebenen-Methode (Demo mit festem Beispiel).

        Beispiel: Maximiere  x₁ + x₂
                  s.t.       2x₁ + 2x₂ ≤ 7
                             x₁, x₂ ≥ 0, ganzzahlig

        LP-Relaxation hat Optimum bei x₁ = x₂ = 1.75 (nicht ganzzahlig).
        Gomory-Schnitt: Füge einen Schnitt hinzu, der x₁+x₂ ≤ 3 erzwingt.

        @return: Dict mit Beschreibung des Verfahrens und Lösung
        @lastModified: 2026-03-10
        """
        # LP-Relaxation des Beispiels lösen: min -(x1+x2) s.t. 2x1+2x2 ≤ 7
        res_relax = linprog(
            c=[-1, -1],
            A_ub=[[2, 2]],
            b_ub=[7],
            bounds=[(0, None), (0, None)],
            method='highs'
        )

        # Gomory-Schnitt manuell: Für x₁+x₂ ≤ floor(7/2)=3
        res_cut = linprog(
            c=[-1, -1],
            A_ub=[[2, 2], [1, 1]],
            b_ub=[7, 3],
            bounds=[(0, None), (0, None)],
            method='highs'
        )

        return {
            'description':      'Gomory-Schnittebenen Demo: max x1+x2 s.t. 2x1+2x2≤7, x ganzzahlig',
            'lp_relaxation_x':  res_relax.x.tolist() if res_relax.success else None,
            'lp_relaxation_val': -res_relax.fun if res_relax.success else None,
            'gomory_cut':       'x1 + x2 ≤ 3',
            'after_cut_x':      res_cut.x.tolist() if res_cut.success else None,
            'after_cut_val':    -res_cut.fun if res_cut.success else None,
            'integer_optimal':  3.0
        }

    @staticmethod
    def knapsack_01(
        values:   List[float],
        weights:  List[float],
        capacity: float
    ) -> Dict:
        """
        0-1-Rucksackproblem (dynamische Programmierung).

        Maximiert den Gesamtwert bei gegebener Kapazitätsbeschränkung:
            max   Σᵢ vᵢ · xᵢ
            s.t.  Σᵢ wᵢ · xᵢ ≤ W
                  xᵢ ∈ {0, 1}

        Zeitkomplexität: O(n · W), Raumkomplexität: O(n · W).

        @param values:   Werte der n Gegenstände
        @param weights:  Gewichte der n Gegenstände
        @param capacity: Maximale Kapazität W (wird zu int gerundet)
        @return: Dict mit 'selected_items', 'total_value', 'total_weight'
        @raises ValueError: Wenn Listen unterschiedliche Längen haben
        @lastModified: 2026-03-10
        """
        n = len(values)
        if len(weights) != n:
            raise ValueError("values und weights müssen gleich lang sein")

        # Kapazität als ganze Zahl (DP-Ansatz erfordert diskrete Kapazität)
        W = int(capacity)

        # DP-Tabelle: dp[i][w] = maximaler Wert mit den ersten i Gegenständen bei Kapazität w
        dp = [[0.0] * (W + 1) for _ in range(n + 1)]

        for i in range(1, n + 1):
            w_i = int(weights[i - 1])
            v_i = values[i - 1]
            for w in range(W + 1):
                # Gegenstand i nicht einpacken
                dp[i][w] = dp[i - 1][w]
                # Gegenstand i einpacken (wenn er passt)
                if w_i <= w:
                    dp[i][w] = max(dp[i][w], dp[i - 1][w - w_i] + v_i)

        # Rückverfolgung: welche Gegenstände wurden ausgewählt?
        selected = []
        w = W
        for i in range(n, 0, -1):
            if dp[i][w] != dp[i - 1][w]:
                selected.append(i - 1)   # 0-basierter Index
                w -= int(weights[i - 1])

        selected.reverse()

        return {
            'selected_items': selected,
            'total_value':    dp[n][W],
            'total_weight':   sum(int(weights[i]) for i in selected)
        }

    @staticmethod
    def knapsack_fractional(
        values:   List[float],
        weights:  List[float],
        capacity: float
    ) -> Dict:
        """
        Fraktionales Rucksackproblem (Greedy-Algorithmus).

        Gegenstände können beliebig geteilt werden.
        Greedy-Strategie: Sortiere nach Wert/Gewicht (absteigend), füge
        ganze Gegenstände ein, bis die Kapazität erschöpft ist.

        Zeitkomplexität: O(n log n) durch Sortierung.
        Optimalität: Greedy liefert immer die optimale Lösung.

        @param values:   Werte der n Gegenstände
        @param weights:  Gewichte der n Gegenstände
        @param capacity: Maximale Kapazität W
        @return: Dict mit 'fractions', 'total_value', 'total_weight'
        @lastModified: 2026-03-10
        """
        n = len(values)
        if len(weights) != n:
            raise ValueError("values und weights müssen gleich lang sein")

        # Wert pro Gewicht berechnen und absteigend sortieren
        items = sorted(
            range(n),
            key=lambda i: values[i] / weights[i] if weights[i] > 0 else math.inf,
            reverse=True
        )

        remaining   = capacity
        total_value = 0.0
        fractions   = [0.0] * n   # Anteil von jedem Gegenstand (0.0 bis 1.0)

        for i in items:
            if remaining <= 0:
                break
            if weights[i] <= remaining:
                # Ganzer Gegenstand passt rein
                fractions[i]  = 1.0
                total_value   += values[i]
                remaining     -= weights[i]
            else:
                # Nur ein Bruchteil passt rein
                frac          = remaining / weights[i]
                fractions[i]  = frac
                total_value   += frac * values[i]
                remaining      = 0.0

        return {
            'fractions':    fractions,
            'total_value':  total_value,
            'total_weight': capacity - remaining
        }

    @staticmethod
    def traveling_salesman_nearest_neighbor(
        dist_matrix: List[List[float]]
    ) -> Dict:
        """
        TSP Nearest-Neighbor-Heuristik.

        Greedy-Heuristik: Beginne in Stadt 0, gehe immer zur nächsten
        noch nicht besuchten Stadt.

        Zeitkomplexität: O(n²).
        Qualität: Im Worst-Case O(log n)-Approximation.

        @param dist_matrix: Entfernungsmatrix (n × n), dist[i][j] = Distanz i→j
        @return: Dict mit 'tour', 'total_distance'
        @raises ValueError: Wenn dist_matrix leer oder nicht quadratisch ist
        @lastModified: 2026-03-10
        """
        D = np.array(dist_matrix, dtype=float)
        n = D.shape[0]
        if n == 0:
            raise ValueError("dist_matrix darf nicht leer sein")
        if D.shape[0] != D.shape[1]:
            raise ValueError("dist_matrix muss quadratisch sein")

        if n == 1:
            return {'tour': [0], 'total_distance': 0.0}

        # Startstadt: 0
        visited    = [False] * n
        tour       = [0]
        visited[0] = True
        total_dist = 0.0

        current = 0
        for _ in range(n - 1):
            # Nächste unbesuchte Stadt finden
            best_dist = math.inf
            best_next = -1
            for j in range(n):
                if not visited[j] and D[current, j] < best_dist:
                    best_dist = D[current, j]
                    best_next = j
            tour.append(best_next)
            visited[best_next] = True
            total_dist += best_dist
            current = best_next

        # Zurück zur Startstadt
        total_dist += D[current, tour[0]]

        return {
            'tour':           tour,
            'total_distance': total_dist
        }

    @staticmethod
    def traveling_salesman_2opt(
        dist_matrix:   List[List[float]],
        initial_tour:  Optional[List[int]] = None
    ) -> Dict:
        """
        2-Opt-Verbesserungsheuristik für das TSP.

        Verbessert eine gegebene Tour durch Vertauschen von 2 Kanten.
        Iteriert, bis keine Verbesserung mehr möglich ist (lokales Optimum).

        Zeitkomplexität pro Iteration: O(n²); Gesamtzahl Iterationen variiert.

        @param dist_matrix:   Entfernungsmatrix (n × n)
        @param initial_tour:  Starttour (Standard: Nearest-Neighbor)
        @return: Dict mit 'tour', 'total_distance', 'iterations'
        @lastModified: 2026-03-10
        """
        D = np.array(dist_matrix, dtype=float)
        n = D.shape[0]
        if n == 0:
            raise ValueError("dist_matrix darf nicht leer sein")

        # Starttour: Nearest-Neighbor wenn nicht angegeben
        if initial_tour is None:
            nn_result  = IntegerProgramming.traveling_salesman_nearest_neighbor(dist_matrix)
            tour       = nn_result['tour']
        else:
            tour = list(initial_tour)

        def tour_length(t: List[int]) -> float:
            """Berechnet die Gesamtlänge einer Tour (inklusive Rückkehr)."""
            length = sum(D[t[i], t[i + 1]] for i in range(len(t) - 1))
            length += D[t[-1], t[0]]
            return length

        improved   = True
        iterations = 0

        while improved:
            improved = False
            iterations += 1
            for i in range(1, n - 1):
                for j in range(i + 1, n):
                    # 2-Opt-Tausch: Segment [i..j] umkehren
                    new_tour = tour[:i] + tour[i:j + 1][::-1] + tour[j + 1:]
                    if tour_length(new_tour) < tour_length(tour) - 1e-10:
                        tour     = new_tour
                        improved = True
                        break
                if improved:
                    break

        return {
            'tour':           tour,
            'total_distance': tour_length(tour),
            'iterations':     iterations
        }


# ===========================================================================
# 3. WARTESCHLANGENTHEORIE
# ===========================================================================

class QueueingTheory:
    """
    Warteschlangentheorie (Queueing Theory).

    Modelle für Ankunfts- und Bedienprozesse:
    - M/M/1: Poisson-Ankünfte, exponentielle Bedienzeiten, 1 Server
    - M/M/s: wie M/M/1, aber s Server
    - M/M/1/K: M/M/1 mit endlicher Kapazität K
    - Erlang-C: Wahrscheinlichkeit, warten zu müssen
    - Erlang-B: Blockierungswahrscheinlichkeit

    Notation:
        λ = Ankunftsrate (Ankünfte pro Zeiteinheit)
        μ = Bedienrate (Bedienungen pro Zeiteinheit pro Server)
        ρ = λ/(s·μ) = Auslastungsgrad

    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """

    @staticmethod
    def mm1_queue(arrival_rate: float, service_rate: float) -> Dict:
        """
        M/M/1-Warteschlange.

        Kenngrößen:
            ρ  = λ/μ            (Auslastung, muss < 1 sein)
            L  = ρ/(1-ρ)        (Mittlere Kundenzahl im System)
            Lq = ρ²/(1-ρ)       (Mittlere Kundenzahl in der Warteschlange)
            W  = 1/(μ-λ)        (Mittlere Verweilzeit im System)
            Wq = λ/(μ(μ-λ))     (Mittlere Wartezeit in der Warteschlange)

        @param arrival_rate:  Ankunftsrate λ (> 0)
        @param service_rate:  Bedienrate μ (> λ für Stabilität)
        @return: Dict mit 'rho', 'L', 'Lq', 'W', 'Wq', 'stable'
        @raises ValueError: Wenn arrival_rate ≤ 0 oder service_rate ≤ 0
        @lastModified: 2026-03-10
        """
        if arrival_rate <= 0:
            raise ValueError("arrival_rate muss > 0 sein")
        if service_rate <= 0:
            raise ValueError("service_rate muss > 0 sein")

        lam = arrival_rate
        mu  = service_rate
        rho = lam / mu

        stable = rho < 1.0

        if not stable:
            # Instabile Schlange: Kenngrößen sind unendlich
            return {
                'rho':    rho,
                'L':      math.inf,
                'Lq':     math.inf,
                'W':      math.inf,
                'Wq':     math.inf,
                'stable': False
            }

        L  = rho / (1 - rho)
        Lq = rho ** 2 / (1 - rho)
        W  = 1.0 / (mu - lam)
        Wq = lam / (mu * (mu - lam))

        return {
            'rho':    rho,
            'L':      L,
            'Lq':     Lq,
            'W':      W,
            'Wq':     Wq,
            'stable': True
        }

    @staticmethod
    def mms_queue(
        arrival_rate:  float,
        service_rate:  float,
        s:             int
    ) -> Dict:
        """
        M/M/s-Warteschlange (s parallele Server).

        Erlang-C-Formel für Wahrscheinlichkeit P₀ (kein Kunde im System):

            P₀ = [Σₙ₌₀ˢ⁻¹ (sρ)ⁿ/n! + (sρ)ˢ/(s!(1-ρ))]⁻¹

        mit ρ = λ/(s·μ) (Auslastung pro Server, muss < 1).

        Kenngrößen:
            Lq = P_C · ρ / (1-ρ)   wobei P_C = Erlang-C-Formel
            L  = Lq + sρ
            Wq = Lq / λ
            W  = Wq + 1/μ

        @param arrival_rate:  Ankunftsrate λ
        @param service_rate:  Bedienrate μ pro Server
        @param s:             Anzahl Server (≥ 1)
        @return: Dict mit 'rho', 'P0', 'Pc', 'L', 'Lq', 'W', 'Wq', 'stable'
        @raises ValueError: Wenn s < 1
        @lastModified: 2026-03-10
        """
        if s < 1:
            raise ValueError("Anzahl Server s muss ≥ 1 sein")
        if arrival_rate <= 0 or service_rate <= 0:
            raise ValueError("arrival_rate und service_rate müssen > 0 sein")

        lam = arrival_rate
        mu  = service_rate
        a   = lam / mu        # Angebotene Verkehrslast (Erlang)
        rho = lam / (s * mu)  # Auslastung pro Server

        stable = rho < 1.0
        if not stable:
            return {
                'rho': rho, 'P0': 0.0, 'Pc': 1.0,
                'L': math.inf, 'Lq': math.inf,
                'W': math.inf, 'Wq': math.inf,
                'stable': False
            }

        # P₀ berechnen
        sum_term = sum((a ** n) / math.factorial(n) for n in range(s))
        last_term = (a ** s) / (math.factorial(s) * (1 - rho))
        P0 = 1.0 / (sum_term + last_term)

        # Erlang-C-Formel: P_C = P(Warten > 0)
        Pc = ((a ** s) / (math.factorial(s) * (1 - rho))) * P0

        # Kenngrößen berechnen
        Lq = Pc * rho / (1 - rho)
        L  = Lq + a         # L = Lq + λ/μ (= Lq + a)
        Wq = Lq / lam
        W  = Wq + 1.0 / mu

        return {
            'rho':    rho,
            'P0':     P0,
            'Pc':     Pc,
            'L':      L,
            'Lq':     Lq,
            'W':      W,
            'Wq':     Wq,
            'stable': True
        }

    @staticmethod
    def mm1k_queue(
        arrival_rate:  float,
        service_rate:  float,
        K:             int
    ) -> Dict:
        """
        M/M/1/K-Warteschlange (endliche Kapazität K).

        Maximal K Kunden im System (einschließlich des Bedienten).
        Ankünfte bei voller Schlange werden abgewiesen.

        Stationäre Verteilung:
            pₙ = (1-ρ) · ρⁿ / (1-ρᴷ⁺¹)  für ρ ≠ 1
            pₙ = 1/(K+1)                  für ρ = 1

        @param arrival_rate:  Ankunftsrate λ
        @param service_rate:  Bedienrate μ
        @param K:             Maximale Systemkapazität (≥ 1)
        @return: Dict mit 'rho', 'P_block', 'L', 'Lq', 'W', 'Wq', 'lambda_eff'
        @raises ValueError: Wenn K < 1
        @lastModified: 2026-03-10
        """
        if K < 1:
            raise ValueError("Systemkapazität K muss ≥ 1 sein")
        if arrival_rate <= 0 or service_rate <= 0:
            raise ValueError("arrival_rate und service_rate müssen > 0 sein")

        lam = arrival_rate
        mu  = service_rate
        rho = lam / mu

        # Stationäre Wahrscheinlichkeiten pₙ
        if abs(rho - 1.0) < 1e-12:
            # Sonderfall ρ = 1
            p = [1.0 / (K + 1)] * (K + 1)
        else:
            # Allgemeiner Fall ρ ≠ 1
            norm = (1 - rho) / (1 - rho ** (K + 1))
            p    = [norm * (rho ** n) for n in range(K + 1)]

        # Blockierungswahrscheinlichkeit: P(System voll) = p_K
        P_block = p[K]

        # Effektive Ankunftsrate (Kunden, die wirklich ankommen)
        lam_eff = lam * (1 - P_block)

        # Mittlere Anzahl im System: L = Σ n·pₙ
        L = sum(n * p[n] for n in range(K + 1))

        # Little's Law: Wq = L/λ_eff - 1/μ
        W  = L / lam_eff if lam_eff > 0 else math.inf
        Wq = W - 1.0 / mu
        Lq = lam_eff * Wq

        return {
            'rho':        rho,
            'P_block':    P_block,
            'L':          L,
            'Lq':         Lq,
            'W':          W,
            'Wq':         Wq,
            'lambda_eff': lam_eff,
            'probs':      p
        }

    @staticmethod
    def erlang_c(
        arrival_rate:  float,
        service_rate:  float,
        s:             int
    ) -> float:
        """
        Erlang-C-Formel: P(Warten > 0) für M/M/s-System.

        Gibt die Wahrscheinlichkeit an, dass ein ankommender Kunde warten muss,
        d.h. alle s Server sind belegt.

            C(s, a) = (aˢ/s!) · (s/(s-a)) / [Σₙ₌₀ˢ⁻¹ aⁿ/n! + (aˢ/s!) · (s/(s-a))]

        mit a = λ/μ (angebotene Last).

        @param arrival_rate:  Ankunftsrate λ
        @param service_rate:  Bedienrate μ pro Server
        @param s:             Anzahl Server
        @return: Erlang-C-Wahrscheinlichkeit ∈ [0, 1]
        @raises ValueError: Wenn s < 1 oder System instabil (ρ ≥ 1)
        @lastModified: 2026-03-10
        """
        if s < 1:
            raise ValueError("s muss ≥ 1 sein")
        result = QueueingTheory.mms_queue(arrival_rate, service_rate, s)
        if not result['stable']:
            return 1.0   # Instables System: immer Wartezeit
        return result['Pc']

    @staticmethod
    def erlang_b(
        arrival_rate:  float,
        service_rate:  float,
        s:             int
    ) -> float:
        """
        Erlang-B-Formel: Blockierungswahrscheinlichkeit für Verlustsystem M/M/s/s.

        Gibt die Wahrscheinlichkeit an, dass alle s Server belegt sind und
        ein ankommender Kunde abgewiesen wird (kein Warten möglich).

        Rekursive Berechnung (numerisch stabil):
            B(0, a) = 1
            B(n, a) = a·B(n-1, a) / (n + a·B(n-1, a))

        @param arrival_rate:  Ankunftsrate λ
        @param service_rate:  Bedienrate μ
        @param s:             Anzahl Server (Kanäle)
        @return: Erlang-B-Wahrscheinlichkeit ∈ [0, 1]
        @raises ValueError: Wenn s < 1
        @lastModified: 2026-03-10
        """
        if s < 1:
            raise ValueError("s muss ≥ 1 sein")
        if arrival_rate <= 0 or service_rate <= 0:
            raise ValueError("Raten müssen > 0 sein")

        a = arrival_rate / service_rate  # Angebotene Last

        # Rekursion für Erlang-B
        B = 1.0
        for n in range(1, s + 1):
            B = (a * B) / (n + a * B)
        return B

    @staticmethod
    def little_law_check(L: float, lambda_: float, W: float) -> Dict:
        """
        Little's Law überprüfen: L = λ · W.

        Little's Law gilt für alle stabilen FIFO-Systeme:
            L = λ · W

        Dabei ist:
            L  = mittlere Anzahl Kunden im System
            λ  = effektive Ankunftsrate
            W  = mittlere Verweilzeit im System

        @param L:        Mittlere Kundenzahl im System
        @param lambda_:  Effektive Ankunftsrate
        @param W:        Mittlere Verweilzeit
        @return: Dict mit 'L_computed', 'L_given', 'relative_error', 'holds'
        @lastModified: 2026-03-10
        """
        L_computed     = lambda_ * W
        relative_error = abs(L_computed - L) / max(abs(L), 1e-10)

        return {
            'L_computed':     L_computed,
            'L_given':        L,
            'relative_error': relative_error,
            'holds':          relative_error < 1e-6
        }


# ===========================================================================
# 4. MARKOV-ENTSCHEIDUNGSPROZESSE
# ===========================================================================

class MarkovDecisionProcess:
    """
    Markov-Entscheidungsprozess (MDP).

    Formale Definition: MDP = (S, A, P, R, γ)
        S  = endliche Zustandsmenge
        A  = endliche Aktionsmenge
        P  = Übergangswahrscheinlichkeiten P(s'|s,a)
        R  = Belohnungsfunktion R(s,a)
        γ  = Diskontierungsfaktor ∈ [0, 1)

    Bellman-Gleichung für optimale Wertfunktion V*:
        V*(s) = max_a [ R(s,a) + γ · Σ_{s'} P(s'|s,a) · V*(s') ]

    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """

    def __init__(
        self,
        states:            List,
        actions:           List,
        transition_probs:  Dict,
        rewards:           Dict,
        gamma:             float = 0.9
    ):
        """
        MDP initialisieren.

        @param states:           Liste der Zustände s ∈ S
        @param actions:          Liste der Aktionen a ∈ A
        @param transition_probs: Dict {(s, a): {s': prob}} – P(s'|s,a)
        @param rewards:          Dict {(s, a): Belohnung} – R(s,a)
        @param gamma:            Diskontierungsfaktor γ ∈ [0, 1)
        @raises ValueError:      Wenn gamma nicht in [0, 1)
        @lastModified: 2026-03-10
        """
        if not (0.0 <= gamma < 1.0):
            raise ValueError(f"gamma muss in [0, 1) liegen, erhalten: {gamma}")

        self.states           = list(states)
        self.actions          = list(actions)
        self.transition_probs = transition_probs
        self.rewards          = rewards
        self.gamma            = gamma

        # Interne Indizes für effiziente Berechnungen
        self.state_index  = {s: i for i, s in enumerate(self.states)}
        self.action_index = {a: i for i, a in enumerate(self.actions)}

        # Wertfunktion und Politik werden bei Bedarf berechnet
        self.V      = {s: 0.0 for s in self.states}
        self.policy = {s: self.actions[0] for s in self.states}

    def bellman_equation(self, V: Dict, s) -> float:
        """
        Bellman-Gleichung für einen Zustand s auswerten.

        Berechnet:
            max_a [ R(s,a) + γ · Σ_{s'} P(s'|s,a) · V(s') ]

        @param V: Aktuelle Wertfunktion {s: Wert}
        @param s: Zustand, für den die Bellman-Gleichung berechnet wird
        @return: Optimaler Bellman-Wert für Zustand s
        @lastModified: 2026-03-10
        """
        best_value = -math.inf

        for a in self.actions:
            # Erwarteter Wert unter Aktion a in Zustand s
            q_value = self.rewards.get((s, a), 0.0)
            transitions = self.transition_probs.get((s, a), {})
            for s_next, prob in transitions.items():
                q_value += self.gamma * prob * V.get(s_next, 0.0)
            best_value = max(best_value, q_value)

        return best_value if best_value > -math.inf else 0.0

    def value_iteration(
        self,
        epsilon:  float = 1e-6,
        max_iter: int   = 1000
    ) -> Dict:
        """
        Wertiteration (Value Iteration).

        Iteriert die Bellman-Optimalitätsgleichung bis zur Konvergenz:
            V_{k+1}(s) = max_a [ R(s,a) + γ · Σ_{s'} P(s'|s,a) · V_k(s') ]

        Konvergenzkriterium: max_s |V_{k+1}(s) - V_k(s)| < ε

        @param epsilon:  Abbruchtoleranz ε (Standard: 1e-6)
        @param max_iter: Maximale Iterationsanzahl
        @return: Dict mit 'value_function', 'policy', 'iterations', 'converged'
        @lastModified: 2026-03-10
        """
        V = {s: 0.0 for s in self.states}

        for iteration in range(max_iter):
            delta = 0.0
            V_new = {}

            for s in self.states:
                v_new      = self.bellman_equation(V, s)
                V_new[s]   = v_new
                delta      = max(delta, abs(v_new - V.get(s, 0.0)))

            V = V_new

            if delta < epsilon:
                # Konvergenz erreicht
                self.V = V
                policy = self.optimal_policy()
                return {
                    'value_function': V,
                    'policy':         policy,
                    'iterations':     iteration + 1,
                    'converged':      True
                }

        # Maximale Iterationsanzahl erreicht
        self.V = V
        policy = self.optimal_policy()
        return {
            'value_function': V,
            'policy':         policy,
            'iterations':     max_iter,
            'converged':      False
        }

    def policy_iteration(self, max_iter: int = 100) -> Dict:
        """
        Politikiteration (Policy Iteration).

        Wechselt zwischen:
        1. Politikbewertung (Policy Evaluation): Berechne V^π
        2. Politikverbesserung (Policy Improvement): Verbessere π greedy

        Konvergiert in endlich vielen Schritten (da |S||A| endlich).

        @param max_iter: Maximale Iterationsanzahl (äußere Schleife)
        @return: Dict mit 'value_function', 'policy', 'iterations', 'converged'
        @lastModified: 2026-03-10
        """
        # Startpolitik: erste Aktion für alle Zustände
        policy = {s: self.actions[0] for s in self.states}
        V      = {s: 0.0 for s in self.states}

        for iteration in range(max_iter):
            # --- Schritt 1: Politikbewertung ---
            # Löse lineares Gleichungssystem: V = R^π + γ·P^π·V
            # Iterative Approximation (einfacher als LGS-Lösung)
            for _ in range(500):
                max_delta = 0.0
                for s in self.states:
                    a         = policy[s]
                    v_new     = self.rewards.get((s, a), 0.0)
                    for s_next, prob in self.transition_probs.get((s, a), {}).items():
                        v_new += self.gamma * prob * V.get(s_next, 0.0)
                    max_delta = max(max_delta, abs(v_new - V[s]))
                    V[s]      = v_new
                if max_delta < 1e-8:
                    break

            # --- Schritt 2: Politikverbesserung ---
            policy_stable = True
            for s in self.states:
                old_action  = policy[s]
                best_action = old_action
                best_value  = -math.inf

                for a in self.actions:
                    q_val       = self.rewards.get((s, a), 0.0)
                    for s_next, prob in self.transition_probs.get((s, a), {}).items():
                        q_val += self.gamma * prob * V.get(s_next, 0.0)
                    if q_val > best_value:
                        best_value  = q_val
                        best_action = a

                policy[s] = best_action
                if best_action != old_action:
                    policy_stable = False

            if policy_stable:
                # Konvergenz: Politik hat sich nicht mehr geändert
                self.V      = V
                self.policy = policy
                return {
                    'value_function': V,
                    'policy':         policy,
                    'iterations':     iteration + 1,
                    'converged':      True
                }

        self.V      = V
        self.policy = policy
        return {
            'value_function': V,
            'policy':         policy,
            'iterations':     max_iter,
            'converged':      False
        }

    def optimal_policy(self) -> Dict:
        """
        Optimale Politik aus der aktuellen Wertfunktion ableiten.

        Greedy-Politik: π*(s) = argmax_a [ R(s,a) + γ · Σ_{s'} P(s'|s,a) · V(s') ]

        @return: Dict {s: optimale Aktion a*}
        @lastModified: 2026-03-10
        """
        policy = {}
        for s in self.states:
            best_action = self.actions[0]
            best_value  = -math.inf

            for a in self.actions:
                q_val = self.rewards.get((s, a), 0.0)
                for s_next, prob in self.transition_probs.get((s, a), {}).items():
                    q_val += self.gamma * prob * self.V.get(s_next, 0.0)
                if q_val > best_value:
                    best_value  = q_val
                    best_action = a

            policy[s] = best_action

        self.policy = policy
        return policy


# ===========================================================================
# 5. NETZWERKFLÜSSE
# ===========================================================================

class NetworkFlows:
    """
    Netzwerkfluss-Algorithmen.

    Implementiert grundlegende Graphalgorithmen für:
    - Max-Flow (Ford-Fulkerson mit BFS = Edmonds-Karp)
    - Min-Cut (über Max-Flow-Min-Cut-Satz)
    - Kürzeste Pfade (Dijkstra, Bellman-Ford)
    - Minimale Spannbäume (Kruskal, Prim)

    Graphrepräsentation:
        Dict {u: {v: Kapazität/Gewicht}} für gerichtete Graphen.

    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """

    @staticmethod
    def max_flow_ford_fulkerson(
        graph:  Dict[int, Dict[int, float]],
        source: int,
        sink:   int
    ) -> Dict:
        """
        Maximaler Fluss via Edmonds-Karp (Ford-Fulkerson + BFS).

        Findet augmentierende Pfade via BFS (Breitensuche) im Residualgraph.
        Zeitkomplexität: O(V · E²)

        Residualgraph: Für jede Kante (u,v) mit Kapazität c und Fluss f:
            - Vorwärtskante (u,v): Restkapazität = c - f
            - Rückwärtskante (v,u): Restkapazität = f

        @param graph:  Adjazenzliste {u: {v: Kapazität}}
        @param source: Quellknoten
        @param sink:   Senkenknoten
        @return: Dict mit 'max_flow', 'flow_dict'
        @lastModified: 2026-03-10
        """
        if not graph:
            return {'max_flow': 0.0, 'flow_dict': {}}

        # Residualgraph aufbauen (Kopie des Graphen mit Rückwärtskanten)
        residual = defaultdict(lambda: defaultdict(float))
        for u in graph:
            for v, cap in graph[u].items():
                residual[u][v] += cap
                # Rückwärtskante mit 0-Kapazität initialisieren (falls nicht vorhanden)
                if v not in residual or u not in residual[v]:
                    residual[v][u] += 0.0

        def bfs_find_path(s: int, t: int) -> Optional[List[int]]:
            """BFS: Augmentierenden Pfad von s nach t im Residualgraph suchen."""
            visited = {s}
            parent  = {s: None}
            queue   = deque([s])

            while queue:
                u = queue.popleft()
                if u == t:
                    # Pfad rekonstruieren
                    path = []
                    while u is not None:
                        path.append(u)
                        u = parent[u]
                    return path[::-1]
                for v in residual[u]:
                    if v not in visited and residual[u][v] > 1e-10:
                        visited.add(v)
                        parent[v] = u
                        queue.append(v)
            return None

        max_flow  = 0.0
        flow_dict = defaultdict(lambda: defaultdict(float))

        # Augmentierende Pfade suchen und Fluss aufaddieren
        while True:
            path = bfs_find_path(source, sink)
            if path is None:
                break  # Kein augmentierender Pfad mehr → Fertig

            # Minimale Restkapazität entlang des Pfades
            path_flow = math.inf
            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]
                path_flow = min(path_flow, residual[u][v])

            # Residualgraph und Fluss aktualisieren
            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]
                residual[u][v] -= path_flow
                residual[v][u] += path_flow
                flow_dict[u][v] += path_flow

            max_flow += path_flow

        return {
            'max_flow':  max_flow,
            'flow_dict': dict(flow_dict)
        }

    @staticmethod
    def min_cut(
        graph:  Dict[int, Dict[int, float]],
        source: int,
        sink:   int
    ) -> Dict:
        """
        Min-Cut via Max-Flow-Min-Cut-Satz.

        Nach dem Max-Flow-Min-Cut-Satz gilt:
            Maximaler Fluss = Kapazität des minimalen Schnitts

        Der minimale Schnitt (S, T) trennt Quelle (∈ S) von Senke (∈ T).
        S = Menge aller Knoten, die im Residualgraph von der Quelle erreichbar sind.

        @param graph:  Adjazenzliste {u: {v: Kapazität}}
        @param source: Quellknoten
        @param sink:   Senkenknoten
        @return: Dict mit 'min_cut_value', 'S_set', 'T_set', 'cut_edges'
        @lastModified: 2026-03-10
        """
        # Maximalen Fluss berechnen (baut Residualgraph intern auf)
        flow_result = NetworkFlows.max_flow_ford_fulkerson(graph, source, sink)

        # Residualgraph nach Max-Flow rekonstruieren
        residual = defaultdict(lambda: defaultdict(float))
        for u in graph:
            for v, cap in graph[u].items():
                residual[u][v] += cap
                residual[v][u] += 0.0

        flow_dict = flow_result['flow_dict']
        for u, neighbors in flow_dict.items():
            for v, f in neighbors.items():
                residual[u][v] -= f
                residual[v][u] += f

        # S = alle von source erreichbaren Knoten im Residualgraph (BFS)
        S_set   = set()
        queue   = deque([source])
        S_set.add(source)
        while queue:
            u = queue.popleft()
            for v in residual[u]:
                if v not in S_set and residual[u][v] > 1e-10:
                    S_set.add(v)
                    queue.append(v)

        # T = alle restlichen Knoten
        all_nodes = set(graph.keys())
        for u in graph:
            all_nodes.update(graph[u].keys())
        T_set = all_nodes - S_set

        # Schnittkanten: (u, v) mit u ∈ S, v ∈ T
        cut_edges = []
        for u in S_set:
            if u in graph:
                for v, cap in graph[u].items():
                    if v in T_set:
                        cut_edges.append((u, v, cap))

        min_cut_value = sum(cap for _, _, cap in cut_edges)

        return {
            'min_cut_value': min_cut_value,
            'S_set':         S_set,
            'T_set':         T_set,
            'cut_edges':     cut_edges
        }

    @staticmethod
    def shortest_path_dijkstra(
        graph:  Dict[int, Dict[int, float]],
        source: int
    ) -> Dict:
        """
        Dijkstra-Algorithmus für kürzeste Pfade.

        Findet kürzeste Pfade von source zu allen anderen Knoten.
        Voraussetzung: Alle Kantengewichte ≥ 0.

        Zeitkomplexität: O((V + E) log V) mit Min-Heap.

        @param graph:  Adjazenzliste {u: {v: Gewicht}} (nur nicht-negative Gewichte)
        @param source: Startknoten
        @return: Dict mit 'distances', 'predecessors'
        @raises ValueError: Wenn negative Gewichte vorhanden sind
        @lastModified: 2026-03-10
        """
        if not graph:
            return {'distances': {source: 0.0}, 'predecessors': {source: None}}

        # Alle Knoten sammeln
        all_nodes = set(graph.keys())
        for u in graph:
            all_nodes.update(graph[u].keys())

        # Distanzen initialisieren
        dist  = {v: math.inf for v in all_nodes}
        dist[source] = 0.0
        pred  = {v: None for v in all_nodes}

        # Min-Heap: (Distanz, Knoten)
        heap = [(0.0, source)]

        while heap:
            d, u = heapq.heappop(heap)

            if d > dist[u]:
                # Veralteter Heap-Eintrag → überspringen
                continue

            for v, w in graph.get(u, {}).items():
                if w < 0:
                    raise ValueError(f"Dijkstra: Negatives Kantengewicht ({u}→{v}: {w})")
                new_dist = dist[u] + w
                if new_dist < dist[v]:
                    dist[v] = new_dist
                    pred[v] = u
                    heapq.heappush(heap, (new_dist, v))

        return {
            'distances':    dist,
            'predecessors': pred
        }

    @staticmethod
    def shortest_path_bellman_ford(
        graph:  Dict[int, Dict[int, float]],
        source: int
    ) -> Dict:
        """
        Bellman-Ford-Algorithmus für kürzeste Pfade (auch mit negativen Gewichten).

        Relaxiert alle Kanten (V-1)-mal. Erkennt negative Zyklen.
        Zeitkomplexität: O(V · E).

        @param graph:  Adjazenzliste {u: {v: Gewicht}}
        @param source: Startknoten
        @return: Dict mit 'distances', 'predecessors', 'negative_cycle'
        @lastModified: 2026-03-10
        """
        if not graph:
            return {'distances': {source: 0.0}, 'predecessors': {source: None}, 'negative_cycle': False}

        # Alle Knoten und Kanten sammeln
        all_nodes = set(graph.keys())
        edges     = []
        for u in graph:
            all_nodes.update(graph[u].keys())
            for v, w in graph[u].items():
                edges.append((u, v, w))

        # Distanzen initialisieren
        dist = {v: math.inf for v in all_nodes}
        dist[source] = 0.0
        pred = {v: None for v in all_nodes}

        # (|V| - 1) Relaxierungsrunden
        n = len(all_nodes)
        for _ in range(n - 1):
            updated = False
            for u, v, w in edges:
                if dist[u] + w < dist[v]:
                    dist[v] = dist[u] + w
                    pred[v] = u
                    updated = True
            if not updated:
                break   # Frühzeitiger Abbruch

        # Negativen Zyklus erkennen: noch eine Runde Relaxierung
        negative_cycle = False
        for u, v, w in edges:
            if dist[u] + w < dist[v] - 1e-10:
                negative_cycle = True
                break

        return {
            'distances':      dist,
            'predecessors':   pred,
            'negative_cycle': negative_cycle
        }

    @staticmethod
    def minimum_spanning_tree_kruskal(
        edges:      List[Tuple[float, int, int]],
        n_vertices: int
    ) -> Dict:
        """
        Kruskal-Algorithmus für minimalen Spannbaum (MST).

        Sortiert Kanten nach Gewicht und fügt gierig hinzu, wenn kein Zyklus entsteht.
        Union-Find-Datenstruktur für effiziente Zykluserkennung.

        Zeitkomplexität: O(E log E).

        @param edges:      Liste von (Gewicht, u, v) – Kanten des ungerichteten Graphen
        @param n_vertices: Anzahl Knoten (0 bis n_vertices-1)
        @return: Dict mit 'mst_edges', 'total_weight', 'n_components'
        @lastModified: 2026-03-10
        """
        if n_vertices == 0:
            return {'mst_edges': [], 'total_weight': 0.0, 'n_components': 0}

        # Union-Find: parent[i] = Repräsentant von i
        parent = list(range(n_vertices))
        rank   = [0] * n_vertices

        def find(x: int) -> int:
            """Pfadkomprimierung: Findet Wurzel von x."""
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x: int, y: int) -> bool:
            """Rang-basierte Vereinigung. Gibt True zurück, wenn x und y verschieden waren."""
            rx, ry = find(x), find(y)
            if rx == ry:
                return False   # Zyklus!
            if rank[rx] < rank[ry]:
                rx, ry = ry, rx
            parent[ry] = rx
            if rank[rx] == rank[ry]:
                rank[rx] += 1
            return True

        # Kanten nach Gewicht sortieren
        sorted_edges = sorted(edges, key=lambda e: e[0])

        mst_edges    = []
        total_weight = 0.0

        for w, u, v in sorted_edges:
            if union(u, v):
                # Kante hinzufügen (kein Zyklus)
                mst_edges.append((w, u, v))
                total_weight += w
                if len(mst_edges) == n_vertices - 1:
                    break   # MST vollständig (n-1 Kanten)

        # Zusammenhangskomponenten zählen
        components = len(set(find(i) for i in range(n_vertices)))

        return {
            'mst_edges':    mst_edges,
            'total_weight': total_weight,
            'n_components': components
        }

    @staticmethod
    def minimum_spanning_tree_prim(adj_matrix: List[List[float]]) -> Dict:
        """
        Prim-Algorithmus für minimalen Spannbaum (MST).

        Wächst den MST iterativ, indem die günstigste Kante zum MST hinzugefügt wird.
        Geeignet für dichte Graphen.

        Zeitkomplexität: O(V²) (ohne Heap) oder O(E log V) (mit Heap).
        Hier: O(V²) Implementierung für Adjazenzmatrizen.

        @param adj_matrix: Adjazenzmatrix (n × n); 0 oder inf = keine Kante
        @return: Dict mit 'mst_edges', 'total_weight'
        @raises ValueError: Wenn adj_matrix leer ist
        @lastModified: 2026-03-10
        """
        A = np.array(adj_matrix, dtype=float)
        n = A.shape[0]
        if n == 0:
            raise ValueError("adj_matrix darf nicht leer sein")

        # 0er-Einträge als fehlende Kanten behandeln (außer Diagonale)
        for i in range(n):
            for j in range(n):
                if i != j and A[i, j] == 0:
                    A[i, j] = math.inf

        in_mst    = [False] * n
        key       = [math.inf] * n   # Minimales Kantengewicht zum MST
        parent    = [-1] * n

        key[0]     = 0.0   # Startknoten: 0
        mst_edges  = []
        total_w    = 0.0

        for _ in range(n):
            # Knoten mit minimalem Key auswählen, der noch nicht im MST ist
            u = -1
            min_key = math.inf
            for v in range(n):
                if not in_mst[v] and key[v] < min_key:
                    min_key = key[v]
                    u = v

            if u == -1:
                break   # Graph nicht zusammenhängend

            in_mst[u] = True

            if parent[u] != -1:
                # Kante zum MST hinzufügen
                mst_edges.append((A[parent[u], u], parent[u], u))
                total_w += A[parent[u], u]

            # Schlüssel der Nachbarn aktualisieren
            for v in range(n):
                if not in_mst[v] and A[u, v] < key[v]:
                    key[v]    = A[u, v]
                    parent[v] = u

        return {
            'mst_edges':    mst_edges,
            'total_weight': total_w
        }


# ===========================================================================
# 6. SPIELTHEORIE
# ===========================================================================

class GameTheory:
    """
    Spieltheorie (Game Theory) für OR-Anwendungen.

    Implementiert grundlegende Konzepte:
    - Nash-Gleichgewicht (Nash Equilibrium)
    - Minimax-Theorem (Nullsummenspiele)
    - Gefangenendilemma
    - Dominante Strategien
    - Pareto-Optimalität

    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """

    @staticmethod
    def nash_equilibrium_2x2(
        A: List[List[float]],
        B: List[List[float]]
    ) -> Dict:
        """
        Nash-Gleichgewicht für 2×2-Bimatrixspiel.

        Findet alle Nash-Gleichgewichte (reine und gemischte Strategien).

        Nash-Gleichgewicht: Strategiepaar (s₁*, s₂*), bei dem kein Spieler
        durch einseitige Abweichung profitiert.

        Für gemischte Strategien (p, 1-p) für Spieler 1 und (q, 1-q) für Spieler 2:
            Spieler 1 indifferent: p so wählen, dass Spieler 2 indifferent ist
                B[0,0]·q + B[1,0]·(1-q) = B[0,1]·q + B[1,1]·(1-q)
            Spieler 2 indifferent: q so wählen, dass Spieler 1 indifferent ist
                A[0,0]·p + A[1,0]·(1-p) = A[0,1]·p + A[1,1]·(1-p)

        @param A: Auszahlungsmatrix für Spieler 1 (2×2)
        @param B: Auszahlungsmatrix für Spieler 2 (2×2)
        @return: Dict mit 'pure_equilibria', 'mixed_equilibrium'
        @lastModified: 2026-03-10
        """
        A = np.array(A, dtype=float)
        B = np.array(B, dtype=float)

        # --- Reine Nash-Gleichgewichte suchen ---
        pure_equilibria = []
        for i in range(2):
            for j in range(2):
                # (i,j) ist NE, wenn:
                # Spieler 1: A[i,j] ≥ A[k,j] für alle k
                # Spieler 2: B[i,j] ≥ B[i,l] für alle l
                p1_best = all(A[i, j] >= A[k, j] for k in range(2))
                p2_best = all(B[i, j] >= B[i, l] for l in range(2))
                if p1_best and p2_best:
                    pure_equilibria.append({
                        'strategies': (i, j),
                        'payoffs':    (float(A[i, j]), float(B[i, j]))
                    })

        # --- Gemischtes Nash-Gleichgewicht ---
        # Spieler 2 mischt: q so dass Spieler 1 indifferent ist
        # A[0,0]·q + A[0,1]·(1-q) = A[1,0]·q + A[1,1]·(1-q)
        # ⟹ q·(A[0,0]-A[0,1]-A[1,0]+A[1,1]) = A[1,1] - A[0,1]
        denom_q = A[0, 0] - A[0, 1] - A[1, 0] + A[1, 1]
        mixed_equilibrium = None

        if abs(denom_q) > 1e-10:
            q = (A[1, 1] - A[0, 1]) / denom_q
            # Spieler 1 mischt: p so dass Spieler 2 indifferent ist
            # B[0,0]·p + B[1,0]·(1-p) = B[0,1]·p + B[1,1]·(1-p)
            # ⟹ p·(B[0,0]-B[0,1]-B[1,0]+B[1,1]) = B[1,1] - B[1,0]
            denom_p = B[0, 0] - B[0, 1] - B[1, 0] + B[1, 1]
            if abs(denom_p) > 1e-10:
                p = (B[1, 1] - B[1, 0]) / denom_p
                if 0 <= p <= 1 and 0 <= q <= 1:
                    # Erwartete Auszahlungen im gemischten Gleichgewicht
                    u1 = p * (A[0, 0] * q + A[0, 1] * (1 - q)) + \
                         (1 - p) * (A[1, 0] * q + A[1, 1] * (1 - q))
                    u2 = q * (B[0, 0] * p + B[1, 0] * (1 - p)) + \
                         (1 - q) * (B[0, 1] * p + B[1, 1] * (1 - p))
                    mixed_equilibrium = {
                        'p':       p,
                        'q':       q,
                        'payoffs': (u1, u2)
                    }

        return {
            'pure_equilibria':  pure_equilibria,
            'mixed_equilibrium': mixed_equilibrium
        }

    @staticmethod
    def minimax_theorem(A: List[List[float]]) -> Dict:
        """
        Minimax-Theorem für Nullsummenspiele (via LP).

        Für jedes endliche Zwei-Personen-Nullsummenspiel gilt:
            max_p min_q pᵀAq = min_q max_p pᵀAq = v* (Spielwert)

        Spieler 1 maximiert (gemischte Strategie p), Spieler 2 minimiert (q).

        LP-Formulierung für Spieler 1 (Maximierung von Spielwert v):
            Substitution yᵢ = pᵢ/(v + shift), v' = 1/(v + shift)
            min  Σᵢ yᵢ   (= 1/v')
            s.t. Aᵀ_shifted · y ≥ 1   (n Bedingungen)
                 yᵢ ≥ 0

        LP-Formulierung für Spieler 2 (Minimierung von Spielwert v):
            Substitution zⱼ = qⱼ/(v + shift)
            max  Σⱼ zⱼ   (= 1/v')  ⟺  min -Σⱼ zⱼ
            s.t. A_shifted · z ≤ 1   (m Bedingungen)
                 zⱼ ≥ 0

        Verschiebung sichert v + shift > 0 (positiver Spielwert).

        @param A: Auszahlungsmatrix für Spieler 1 (m × n)
        @return: Dict mit 'game_value', 'p_strategy', 'q_strategy', 'success'
        @lastModified: 2026-03-10
        """
        A_np = np.array(A, dtype=float)
        m, n = A_np.shape

        # Verschiebung: Spielwert nach shift > 0 (garantiert positive Schranke)
        shift     = abs(A_np.min()) + 1.0
        A_shifted = A_np + shift

        # --- LP für Spieler 1: min Σyᵢ  s.t. Aᵀ_shifted · y ≥ 1, y ≥ 0 ---
        # ⟺  min Σyᵢ  s.t. -Aᵀ_shifted · y ≤ -1, y ≥ 0
        c_p     = [1.0] * m
        A_ub_p  = (-A_shifted.T).tolist()   # (n × m)-Matrix
        b_ub_p  = [-1.0] * n
        bounds_p = [(0, None)] * m

        res_p = linprog(c=c_p, A_ub=A_ub_p, b_ub=b_ub_p, bounds=bounds_p, method='highs')

        if not res_p.success:
            return {'game_value': None, 'p_strategy': None, 'q_strategy': None, 'success': False}

        # Strategie p aus y rückgewinnen
        y_p     = res_p.x
        sum_y_p = y_p.sum()
        p_opt   = (y_p / sum_y_p).tolist() if sum_y_p > 1e-10 else (y_p / max(y_p.sum(), 1e-10)).tolist()
        # Spielwert aus Summe: v = 1/Σyᵢ - shift
        v_star  = (1.0 / sum_y_p) - shift if sum_y_p > 1e-10 else None

        # --- LP für Spieler 2: min -Σzⱼ  s.t. A_shifted · z ≤ 1, z ≥ 0 ---
        c_q      = [-1.0] * n
        A_ub_q   = A_shifted.tolist()   # (m × n)-Matrix
        b_ub_q   = [1.0] * m
        bounds_q = [(0, None)] * n

        res_q = linprog(c=c_q, A_ub=A_ub_q, b_ub=b_ub_q, bounds=bounds_q, method='highs')

        q_opt = None
        if res_q.success:
            z_q     = res_q.x
            sum_z_q = z_q.sum()
            q_opt   = (z_q / sum_z_q).tolist() if sum_z_q > 1e-10 else z_q.tolist()

        return {
            'game_value': float(v_star) if v_star is not None else None,
            'p_strategy': p_opt,
            'q_strategy': q_opt,
            'success':    True
        }

    @staticmethod
    def prisoners_dilemma() -> Dict:
        """
        Gefangenendilemma-Analyse.

        Klassisches Beispiel für ein Nicht-Nullsummenspiel:
        Beide Spieler können kooperieren (C) oder defektieren (D).

        Auszahlungsmatrix (Jahre Gefängnis, negiert als Nutzen):
                    C        D
            C  (-1, -1)  (-3,  0)
            D  ( 0, -3)  (-2, -2)

        Dominante Strategie: D (defektieren) für beide Spieler.
        Nash-Gleichgewicht: (D, D) mit Auszahlung (-2, -2).
        Pareto-Optimum: (C, C) mit Auszahlung (-1, -1) → besser für beide!

        @return: Dict mit Spielmatrizen, NE, Pareto-Optimum und Analyse
        @lastModified: 2026-03-10
        """
        # Auszahlungsmatrizen (negierte Jahre = Nutzen)
        A = [[-1, -3], [0, -2]]   # Spieler 1
        B = [[-1, 0], [-3, -2]]   # Spieler 2

        ne_result = GameTheory.nash_equilibrium_2x2(A, B)

        return {
            'payoff_matrix_A':    A,
            'payoff_matrix_B':    B,
            'strategy_labels':    ['C (kooperieren)', 'D (defektieren)'],
            'nash_equilibrium':   ne_result['pure_equilibria'],
            'dominant_strategy':  'D (defektieren) für beide Spieler',
            'pareto_optimum':     {'strategies': (0, 0), 'payoffs': (-1, -1)},
            'dilemma_explanation': (
                'Obwohl (C,C) für beide besser ist, ist (D,D) das NE, '
                'weil D eine dominante Strategie für jeden Spieler ist.'
            )
        }

    @staticmethod
    def dominant_strategy(payoff_matrix: List[List[float]]) -> Dict:
        """
        Dominante Strategie finden.

        Eine Strategie i dominiert Strategie k (strikt), wenn:
            Aᵢⱼ > Aₖⱼ für alle j

        Eine Strategie i dominiert Strategie k (schwach), wenn:
            Aᵢⱼ ≥ Aₖⱼ für alle j (mit mindestens einer strikten Ungleichung)

        @param payoff_matrix: Auszahlungsmatrix (m × n); Zeilen = Strategien des Spielers
        @return: Dict mit 'strictly_dominant', 'weakly_dominant', 'dominated_strategies'
        @lastModified: 2026-03-10
        """
        P = np.array(payoff_matrix, dtype=float)
        m = P.shape[0]

        strictly_dominant  = []   # Strategien, die alle anderen strikt dominieren
        weakly_dominant    = []   # Strategien, die alle anderen schwach dominieren
        dominated          = []   # Dominierte Strategien

        for i in range(m):
            is_strictly = True
            is_weakly   = True
            is_dominated = False

            for k in range(m):
                if k == i:
                    continue
                # Prüfe ob i k strikt dominiert
                if not np.all(P[i] > P[k]):
                    is_strictly = False
                # Prüfe ob i k schwach dominiert
                if not (np.all(P[i] >= P[k]) and np.any(P[i] > P[k])):
                    is_weakly = False
                # Prüfe ob i durch k dominiert wird
                if np.all(P[k] >= P[i]) and np.any(P[k] > P[i]):
                    is_dominated = True

            if is_strictly:
                strictly_dominant.append(i)
            if is_weakly:
                weakly_dominant.append(i)
            if is_dominated:
                dominated.append(i)

        return {
            'strictly_dominant':  strictly_dominant,
            'weakly_dominant':    weakly_dominant,
            'dominated_strategies': dominated
        }

    @staticmethod
    def pareto_optimal(outcomes: List[Tuple[float, ...]]) -> Dict:
        """
        Pareto-optimale Auszahlungen bestimmen.

        Ein Ergebnis x ist Pareto-optimal (Pareto-effizient), wenn es kein
        anderes Ergebnis y gibt, das mindestens einen Spieler besser und
        keinen schlechter stellt.

        Formal: x ist Pareto-optimal, wenn kein y mit:
            yᵢ ≥ xᵢ für alle i  und  yⱼ > xⱼ für mindestens ein j

        @param outcomes: Liste von Auszahlungstupeln (u₁, u₂, ..., uₖ) pro Spieler
        @return: Dict mit 'pareto_optimal', 'dominated_outcomes'
        @lastModified: 2026-03-10
        """
        n = len(outcomes)
        is_pareto = [True] * n

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                # Prüfe ob j outcome i dominiert (j besser in allen oder gleich gut, besser in min. einem)
                if all(outcomes[j][k] >= outcomes[i][k] for k in range(len(outcomes[i]))) and \
                   any(outcomes[j][k] > outcomes[i][k] for k in range(len(outcomes[i]))):
                    is_pareto[i] = False
                    break

        pareto_optimal_outcomes = [outcomes[i] for i in range(n) if is_pareto[i]]
        dominated_outcomes      = [outcomes[i] for i in range(n) if not is_pareto[i]]

        return {
            'pareto_optimal':    pareto_optimal_outcomes,
            'dominated_outcomes': dominated_outcomes,
            'pareto_indices':    [i for i in range(n) if is_pareto[i]]
        }


# ===========================================================================
# 7. STOCHASTISCHE OPTIMIERUNG
# ===========================================================================

class StochasticOptimization:
    """
    Stochastische Optimierung.

    Methoden für Optimierungsprobleme unter Unsicherheit:
    - Stochastischer Gradientenabstieg (SGD)
    - Simulated Annealing (SA)
    - Zeitungshändlerproblem (Newsvendor)
    - Stochastische Programmierung (Demo)

    @author: Kurt Ingwer
    @lastModified: 2026-03-10
    """

    @staticmethod
    def stochastic_gradient_descent(
        f:        Callable[[np.ndarray], float],
        grad_f:   Callable[[np.ndarray, int], np.ndarray],
        x0:       List[float],
        n_steps:  int,
        lr:       float = 0.01,
        seed:     int   = 42
    ) -> Dict:
        """
        Stochastischer Gradientenabstieg (SGD).

        Aktualisierungsregel:
            x_{t+1} = x_t - α · ∇f_i(x_t)

        wobei i zufällig aus {0, ..., n_steps-1} gewählt wird.

        @param f:        Zielfunktion f(x) → float
        @param grad_f:   Stochastischer Gradient ∇f_i(x, i) → np.ndarray
        @param x0:       Startpunkt
        @param n_steps:  Anzahl SGD-Schritte
        @param lr:       Lernrate α (Standard: 0.01)
        @param seed:     Zufallsseed für Reproduzierbarkeit
        @return: Dict mit 'x_opt', 'f_opt', 'history'
        @lastModified: 2026-03-10
        """
        rng = np.random.default_rng(seed)
        x   = np.array(x0, dtype=float)

        history = [float(f(x))]

        for t in range(n_steps):
            # Zufälligen Datenindex wählen
            i  = rng.integers(0, n_steps)
            g  = grad_f(x, int(i))
            # SGD-Update
            x -= lr * g
            if t % max(1, n_steps // 20) == 0:
                history.append(float(f(x)))

        return {
            'x_opt':   x.tolist(),
            'f_opt':   float(f(x)),
            'history': history
        }

    @staticmethod
    def simulated_annealing(
        energy_fn: Callable[[np.ndarray], float],
        x0:        List[float],
        T_init:    float = 1.0,
        cooling:   float = 0.95,
        n_steps:   int   = 1000,
        seed:      int   = 42
    ) -> Dict:
        """
        Simulated Annealing (Simuliertes Ausglühen).

        Metaheuristik zur globalen Optimierung. Akzeptiert schlechtere Lösungen
        mit Wahrscheinlichkeit exp(-ΔE / T) (Boltzmann-Akzeptanz).

        Temperaturabkühlungsplan:
            T(k) = T₀ · cooling^k

        @param energy_fn: Energiefunktion E(x) → float (wird minimiert)
        @param x0:        Startlösung
        @param T_init:    Anfangstemperatur T₀
        @param cooling:   Abkühlungsrate ∈ (0, 1)
        @param n_steps:   Anzahl Simulationsschritte
        @param seed:      Zufallsseed
        @return: Dict mit 'x_opt', 'energy_opt', 'T_final', 'history'
        @lastModified: 2026-03-10
        """
        rng      = np.random.default_rng(seed)
        x        = np.array(x0, dtype=float)
        E        = energy_fn(x)
        x_best   = x.copy()
        E_best   = E
        T        = T_init
        history  = [E]

        for step in range(n_steps):
            # Nachbar erzeugen: kleiner zufälliger Schritt
            x_new = x + rng.normal(0, T, size=x.shape)
            E_new = energy_fn(x_new)

            delta_E = E_new - E

            # Metropolis-Akzeptanzkriterium
            if delta_E < 0 or rng.random() < math.exp(-delta_E / max(T, 1e-10)):
                x = x_new
                E = E_new

            # Beste Lösung verfolgen
            if E < E_best:
                E_best = E
                x_best = x.copy()

            # Temperatur abkühlen
            T *= cooling

            if step % max(1, n_steps // 20) == 0:
                history.append(E_best)

        return {
            'x_opt':       x_best.tolist(),
            'energy_opt':  E_best,
            'T_final':     T,
            'history':     history
        }

    @staticmethod
    def newsvendor_problem(
        demand_samples:  List[float],
        cost_overstock:  float,
        cost_understock: float
    ) -> Dict:
        """
        Zeitungshändlerproblem (Newsvendor Problem).

        Optimale Bestellmenge Q* minimiert erwartete Kosten:
            E[C(Q)] = c_o · E[max(Q-D, 0)] + c_u · E[max(D-Q, 0)]

        Optimale Lösung via kritischer Ratio (CR):
            Q* = F⁻¹(CR)  mit  CR = c_u / (c_u + c_o)

        F⁻¹ = empirisches Quantil der Nachfrageverteilung.

        @param demand_samples:  Historische Nachfragedaten
        @param cost_overstock:  Kosten pro Einheit bei Überbestand (c_o > 0)
        @param cost_understock: Kosten pro Einheit bei Unterbestand (c_u > 0)
        @return: Dict mit 'optimal_quantity', 'critical_ratio', 'expected_cost'
        @raises ValueError: Wenn cost_overstock oder cost_understock ≤ 0
        @lastModified: 2026-03-10
        """
        if cost_overstock <= 0 or cost_understock <= 0:
            raise ValueError("Kosten müssen > 0 sein")

        c_o = cost_overstock
        c_u = cost_understock
        D   = np.array(demand_samples, dtype=float)

        # Kritische Ratio
        CR = c_u / (c_u + c_o)

        # Optimale Bestellmenge: empirisches CR-Quantil
        Q_star = float(np.quantile(D, CR))

        # Erwartete Kosten bei Q*
        overstock_cost  = c_o * np.mean(np.maximum(Q_star - D, 0))
        understock_cost = c_u * np.mean(np.maximum(D - Q_star, 0))
        expected_cost   = overstock_cost + understock_cost

        return {
            'optimal_quantity': Q_star,
            'critical_ratio':   CR,
            'expected_cost':    expected_cost,
            'overstock_cost':   overstock_cost,
            'understock_cost':  understock_cost
        }

    @staticmethod
    def stochastic_programming_demo(
        scenarios:     List[Dict],
        probabilities: List[float]
    ) -> Dict:
        """
        Stochastische Programmierung: Zwei-Stufen-Modell (Demo).

        Erwartungswert-Modell:
            min  cᵀx + Σₛ pₛ · qₛᵀyₛ
            s.t. Ax = b               (Erste-Stufe-Constraint)
                 Tₛx + Wyₛ = hₛ       (Zweite-Stufe-Constraint je Szenario s)
                 x ≥ 0, yₛ ≥ 0

        Hier vereinfacht: Minimiere Erwartungswert der Zielfunktion über Szenarien.
        Szenarien haben Format: {'c': Kosten, 'A_ub': Matrix, 'b_ub': Vektor}

        @param scenarios:     Liste von Szenario-Dicts mit LP-Parametern
        @param probabilities: Wahrscheinlichkeiten der Szenarien (Summe = 1)
        @return: Dict mit 'expected_value', 'scenario_results', 'wait_and_see_value'
        @lastModified: 2026-03-10
        """
        if abs(sum(probabilities) - 1.0) > 1e-6:
            raise ValueError("Summe der Wahrscheinlichkeiten muss 1 sein")

        scenario_results = []
        wait_and_see     = 0.0   # Perfekte Information

        for scenario, prob in zip(scenarios, probabilities):
            c     = scenario['c']
            A_ub  = scenario.get('A_ub', None)
            b_ub  = scenario.get('b_ub', None)
            n     = len(c)
            bounds = [(0, None)] * n

            res = linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')
            if res.success:
                scenario_results.append({
                    'probability':     prob,
                    'optimal_value':   res.fun,
                    'optimal_x':       res.x.tolist()
                })
                wait_and_see += prob * res.fun
            else:
                scenario_results.append({
                    'probability':   prob,
                    'optimal_value': None,
                    'message':       res.message
                })

        # Erwartungswert bei bekannter Strategie
        expected_value = sum(
            r['probability'] * r['optimal_value']
            for r in scenario_results
            if r['optimal_value'] is not None
        )

        return {
            'expected_value':       expected_value,
            'wait_and_see_value':   wait_and_see,
            'evpi':                 abs(expected_value - wait_and_see),
            'scenario_results':     scenario_results
        }
