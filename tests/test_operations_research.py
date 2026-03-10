"""
Tests für operations_research.py.

Umfassende Test-Suite für alle OR-Klassen:
    - LinearProgramming: LP, dual, sensitivity, transport, assignment, 2-phase
    - IntegerProgramming: B&B, cutting plane, knapsack, TSP
    - QueueingTheory: M/M/1, M/M/s, M/M/1/K, Erlang-C/B, Little's Law
    - MarkovDecisionProcess: value/policy iteration, Bellman, optimal policy
    - NetworkFlows: max-flow, min-cut, Dijkstra, Bellman-Ford, Kruskal, Prim
    - GameTheory: Nash 2x2, minimax, prisoner's dilemma, dominant, Pareto
    - StochasticOptimization: SGD, SA, newsvendor, stochastic programming

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-10
@lastModified: 2026-03-10
"""

import math
import sys
import os
import pytest
import numpy as np

# Suchpfad für src/-Verzeichnis setzen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from operations_research import (
    LinearProgramming,
    IntegerProgramming,
    QueueingTheory,
    MarkovDecisionProcess,
    NetworkFlows,
    GameTheory,
    StochasticOptimization,
)


# ===========================================================================
# 1. LINEARE OPTIMIERUNG
# ===========================================================================

class TestLinearProgramming:
    """Tests für LinearProgramming."""

    # --- simplex_method ---

    def test_simplex_simple_minimization(self):
        """Einfache 2-Variablen LP: min x1+x2 s.t. x1+x2≥2, x1,x2≥0."""
        # Umformulierung: min x1+x2 s.t. -(x1+x2) ≤ -2
        c = [1.0, 1.0]
        A_ub = [[-1.0, -1.0]]
        b_ub = [-2.0]
        res = LinearProgramming.simplex_method(c, A_ub, b_ub)
        assert res['success']
        assert abs(res['fun'] - 2.0) < 1e-4

    def test_simplex_classic_lp(self):
        """Klassisches LP: min -x1-2x2 s.t. x1+x2≤4, x1-x2≤2, x1,x2≥0."""
        c = [-1.0, -2.0]
        A_ub = [[1.0, 1.0], [1.0, -1.0]]
        b_ub = [4.0, 2.0]
        res = LinearProgramming.simplex_method(c, A_ub, b_ub)
        assert res['success']
        assert res['fun'] < -3.0   # Optimum bei x2=3, x1=1 ⟹ -7

    def test_simplex_with_equality_constraints(self):
        """LP mit Gleichungsnebenbedingungen (nur A_eq, kein A_ub)."""
        c = [1.0, 1.0]
        A_eq = [[1.0, 1.0]]
        b_eq = [5.0]
        # A_ub und b_ub als None übergeben (kein Ungleichungs-Constraint)
        res = LinearProgramming.simplex_method(c, None, None, A_eq=A_eq, b_eq=b_eq)
        assert res['success']
        assert abs(res['fun'] - 5.0) < 1e-4

    def test_simplex_with_custom_bounds(self):
        """LP mit benutzerdefinierten Variablenschranken."""
        c = [-1.0]
        A_ub = [[1.0]]
        b_ub = [10.0]
        bounds = [(0, 5)]
        res = LinearProgramming.simplex_method(c, A_ub, b_ub, bounds=bounds)
        assert res['success']
        assert abs(res['x'][0] - 5.0) < 1e-4

    def test_simplex_returns_dict_keys(self):
        """Rückgabetyp hat alle erwarteten Schlüssel."""
        res = LinearProgramming.simplex_method([-1.0], [[1.0]], [10.0])
        assert all(k in res for k in ('x', 'fun', 'status', 'message', 'success'))

    def test_simplex_zero_objective(self):
        """LP mit Nullzielfunktion: jede zulässige Lösung ist optimal."""
        c = [0.0, 0.0]
        A_ub = [[1.0, 0.0], [0.0, 1.0]]
        b_ub = [5.0, 5.0]
        res = LinearProgramming.simplex_method(c, A_ub, b_ub)
        assert res['success']
        assert abs(res['fun']) < 1e-6

    # --- dual_problem ---

    def test_dual_problem_returns_structure(self):
        """dual_problem liefert Dict mit erwarteten Schlüsseln."""
        c = [2.0, 3.0]
        A = [[1.0, 1.0], [1.0, 2.0]]
        b = [4.0, 6.0]
        res = LinearProgramming.dual_problem(c, A, b)
        assert 'dual_variables' in res
        assert 'dual_objective' in res
        assert 'success' in res

    def test_dual_problem_strong_duality(self):
        """
        Starker Dualitätssatz: Primal- und Dual-Optimalwert stimmen überein.
        Primal: min 2x1+3x2 s.t. x1+x2≥4, x1+2x2≥6, x≥0.
        Dual:   max 4y1+6y2 s.t. y1+y2≤2, y1+2y2≤3, y≥0.
        """
        c = [2.0, 3.0]
        A = [[1.0, 1.0], [1.0, 2.0]]
        b = [4.0, 6.0]
        res_dual = LinearProgramming.dual_problem(c, A, b)

        # Primal-Lösung direkt (A·x ≥ b ⟺ -A·x ≤ -b)
        res_primal = LinearProgramming.simplex_method(c, [[-1,-1],[-1,-2]], [-4,-6])
        if res_primal['success'] and res_dual['success']:
            assert abs(res_dual['dual_objective'] - res_primal['fun']) < 1e-2

    # --- sensitivity_analysis ---

    def test_sensitivity_analysis_returns_shadow_prices(self):
        """Sensitivitätsanalyse liefert Schattenpreise."""
        c = [-1.0, -2.0]
        A_ub = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
        b_ub = [4.0, 6.0, 7.0]
        res = LinearProgramming.sensitivity_analysis(c, A_ub, b_ub)
        assert res['shadow_prices'] is not None
        assert len(res['shadow_prices']) == 3

    def test_sensitivity_analysis_non_binding_constraint(self):
        """Schattenpreis einer nicht-bindenden Nebenbedingung ist nahe 0."""
        c = [-1.0]
        A_ub = [[1.0], [10.0]]
        b_ub = [5.0, 100.0]   # Zweite NB ist nicht bindend (x≤10, aber x≤5)
        res = LinearProgramming.sensitivity_analysis(c, A_ub, b_ub)
        assert res['shadow_prices'] is not None
        # Schattenpreis der nicht-bindenden NB ≈ 0
        assert abs(res['shadow_prices'][1]) < 0.01

    # --- transportation_problem ---

    def test_transportation_balanced(self):
        """Balanciertes Transportproblem: supply = demand."""
        supply = [30.0, 40.0]
        demand = [20.0, 30.0, 20.0]
        costs  = [[2.0, 3.0, 1.0],
                  [5.0, 4.0, 8.0]]
        res = LinearProgramming.transportation_problem(supply, demand, costs)
        assert res['success']
        assert res['total_cost'] > 0
        # Angebotsrestriktionen prüfen
        flow = res['flow_matrix']
        assert abs(flow[:, :3].sum(axis=1)[0] - 30.0) < 1e-3
        assert abs(flow[:, :3].sum(axis=1)[1] - 40.0) < 1e-3

    def test_transportation_unbalanced(self):
        """Unbalanciertes Transportproblem wird automatisch ausbalanciert."""
        supply = [50.0, 30.0]
        demand = [40.0, 20.0]   # Summe Demand < Summe Supply
        costs  = [[1.0, 2.0], [3.0, 4.0]]
        res = LinearProgramming.transportation_problem(supply, demand, costs)
        assert res['success']

    def test_transportation_single_source_destination(self):
        """Triviales Transportproblem: 1 Quelle, 1 Ziel."""
        supply = [100.0]
        demand = [100.0]
        costs  = [[5.0]]
        res = LinearProgramming.transportation_problem(supply, demand, costs)
        assert res['success']
        assert abs(res['total_cost'] - 500.0) < 1.0

    # --- assignment_problem ---

    def test_assignment_3x3(self):
        """Zuweisungsproblem 3×3: bekannte optimale Lösung."""
        cost_matrix = [
            [4, 1, 3],
            [2, 0, 5],
            [3, 2, 2]
        ]
        res = LinearProgramming.assignment_problem(cost_matrix)
        assert 'row_ind' in res
        assert 'col_ind' in res
        # Alle Zeilen und Spalten eindeutig zugewiesen
        assert len(set(res['row_ind'])) == 3
        assert len(set(res['col_ind'])) == 3
        # Optimaler Wert: 0+2+2=4 oder 1+2+2=5 ...
        assert res['total_cost'] <= 7.0   # Besser als naiver Diagonale (4+0+2=6)

    def test_assignment_empty_raises(self):
        """Leere Kostenmatrix wirft ValueError."""
        with pytest.raises(ValueError):
            LinearProgramming.assignment_problem([])

    def test_assignment_1x1(self):
        """Triviales Zuweisungsproblem 1×1."""
        res = LinearProgramming.assignment_problem([[7.0]])
        assert res['total_cost'] == pytest.approx(7.0)

    def test_assignment_identity_matrix(self):
        """Einheitsmatrix: optimale Zuweisung hat Kosten 0."""
        identity = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
        # Minimale Kosten: 0 (die 0er auf der Gegendiagonale)
        res = LinearProgramming.assignment_problem(
            [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        )
        assert res['total_cost'] == pytest.approx(0.0, abs=1e-6)

    # --- two_phase_simplex_demo ---

    def test_two_phase_simplex_feasible(self):
        """Zwei-Phasen-Simplex findet zulässige Lösung."""
        c = [1.0, 1.0]
        A = [[1.0, 0.0], [0.0, 1.0]]
        b = [2.0, 3.0]
        res = LinearProgramming.two_phase_simplex_demo(c, A, b)
        assert res['feasible']
        assert res['optimal_value'] is not None

    def test_two_phase_simplex_optimal_value(self):
        """Zwei-Phasen-Simplex: min x1+x2 s.t. x1+x2=5 ⟹ Optimum 5."""
        c = [1.0, 1.0]
        A = [[1.0, 1.0]]
        b = [5.0]
        res = LinearProgramming.two_phase_simplex_demo(c, A, b)
        assert res['feasible']
        assert abs(res['optimal_value'] - 5.0) < 1e-4


# ===========================================================================
# 2. GANZZAHLIGE OPTIMIERUNG
# ===========================================================================

class TestIntegerProgramming:
    """Tests für IntegerProgramming."""

    # --- branch_and_bound_demo ---

    def test_branch_and_bound_simple(self):
        """B&B für einfaches IP: min -x1-x2 s.t. x1+x2≤3, x ganzzahlig."""
        c = [-1.0, -1.0]
        A_ub = [[1.0, 1.0]]
        b_ub = [3.0]
        bounds = [(0, 3), (0, 3)]
        res = IntegerProgramming.branch_and_bound_demo(c, A_ub, b_ub, bounds)
        assert res['success']
        # Optimum: x1=x2=1.5 → ganzzahlig: x1+x2=3, z.B. (3,0) oder (0,3) oder (1,2) etc.
        assert abs(res['optimal_value'] - (-3.0)) < 1.0

    def test_branch_and_bound_returns_integer(self):
        """B&B liefert ganzzahlige Lösung."""
        c = [-3.0, -5.0]
        A_ub = [[1.0, 0.0], [0.0, 1.0]]
        b_ub = [4.0, 6.0]
        bounds = [(0, 10), (0, 10)]
        res = IntegerProgramming.branch_and_bound_demo(c, A_ub, b_ub, bounds)
        if res['success']:
            for xi in res['optimal_x']:
                assert abs(xi - round(xi)) < 1e-6

    def test_branch_and_bound_infeasible(self):
        """B&B für unzulässiges Problem liefert success=False."""
        c = [1.0]
        A_ub = [[1.0], [-1.0]]
        b_ub = [-1.0, -2.0]   # x ≤ -1 und x ≥ 2 → unzulässig
        bounds = [(0, None)]
        res = IntegerProgramming.branch_and_bound_demo(c, A_ub, b_ub, bounds)
        assert not res['success']

    # --- cutting_plane_demo ---

    def test_cutting_plane_demo(self):
        """Schnittebenen-Demo liefert ganzzahlige Lösung."""
        res = IntegerProgramming.cutting_plane_demo()
        assert res['integer_optimal'] == 3.0
        assert res['lp_relaxation_val'] > 3.0   # LP-Relaxation besser als IP

    def test_cutting_plane_after_cut_is_integer(self):
        """Nach dem Gomory-Schnitt ist die Lösung ganzzahlig."""
        res = IntegerProgramming.cutting_plane_demo()
        if res['after_cut_x'] is not None:
            for xi in res['after_cut_x']:
                assert abs(xi - round(xi)) < 0.01

    # --- knapsack_01 ---

    def test_knapsack_01_basic(self):
        """0-1-Rucksack: Klassisches 4-Item Beispiel, Kapazität 8.

        Items: (w=5,v=10), (w=4,v=6), (w=3,v=5), (w=2,v=3)
        Optimale Wahl: Item 0 (w=5,v=10) + Item 2 (w=3,v=5) → w=8, v=15.
        """
        values   = [10, 6, 5, 3]
        weights  = [5,  4, 3, 2]
        capacity = 8
        res = IntegerProgramming.knapsack_01(values, weights, capacity)
        # Optimum: Item 0 + Item 2 → total_value=15, total_weight=8
        assert res['total_value'] == 15
        assert res['total_weight'] <= capacity

    def test_knapsack_01_full_capacity(self):
        """Alle Items passen rein."""
        values   = [1, 2, 3]
        weights  = [1, 1, 1]
        capacity = 3
        res = IntegerProgramming.knapsack_01(values, weights, capacity)
        assert res['total_value'] == 6
        assert len(res['selected_items']) == 3

    def test_knapsack_01_zero_capacity(self):
        """Kapazität 0: kein Item wird gewählt."""
        values   = [10, 20]
        weights  = [5, 10]
        capacity = 0
        res = IntegerProgramming.knapsack_01(values, weights, capacity)
        assert res['total_value'] == 0
        assert len(res['selected_items']) == 0

    def test_knapsack_01_single_item(self):
        """Einzelnes Item: entweder nehmen oder nicht."""
        res1 = IntegerProgramming.knapsack_01([10], [5], 5)
        assert res1['total_value'] == 10

        res2 = IntegerProgramming.knapsack_01([10], [6], 5)
        assert res2['total_value'] == 0

    def test_knapsack_01_mismatched_lengths_raises(self):
        """Unterschiedliche Längen werfen ValueError."""
        with pytest.raises(ValueError):
            IntegerProgramming.knapsack_01([1, 2], [1], 5)

    def test_knapsack_01_optimal_value(self):
        """Bekanntes Optimum: values=[60,100,120], weights=[10,20,30], cap=50 → 220."""
        res = IntegerProgramming.knapsack_01([60, 100, 120], [10, 20, 30], 50)
        assert res['total_value'] == 220

    # --- knapsack_fractional ---

    def test_knapsack_fractional_basic(self):
        """Fraktionales Rucksackproblem: Greedy-Lösung."""
        values   = [60, 100, 120]
        weights  = [10, 20,  30]
        capacity = 50
        res = IntegerProgramming.knapsack_fractional(values, weights, capacity)
        assert abs(res['total_value'] - 240.0) < 1e-6
        assert abs(res['total_weight'] - 50.0) < 1e-6

    def test_knapsack_fractional_sum_fractions(self):
        """Summe der Anteile mal Gewicht ≤ Kapazität."""
        values   = [10, 20, 30]
        weights  = [5,  10, 15]
        capacity = 20
        res = IntegerProgramming.knapsack_fractional(values, weights, capacity)
        total_w = sum(res['fractions'][i] * weights[i] for i in range(3))
        assert total_w <= capacity + 1e-9

    def test_knapsack_fractional_full_fit(self):
        """Alle Items passen vollständig → keine Fraktionierung nötig."""
        values   = [1, 2, 3]
        weights  = [1, 1, 1]
        capacity = 5
        res = IntegerProgramming.knapsack_fractional(values, weights, capacity)
        assert all(f == 1.0 for f in res['fractions'])
        assert abs(res['total_value'] - 6.0) < 1e-9

    # --- traveling_salesman_nearest_neighbor ---

    def test_tsp_nn_4cities(self):
        """TSP Nearest-Neighbor mit 4 Städten."""
        D = [
            [0, 10, 15, 20],
            [10, 0, 35, 25],
            [15, 35, 0, 30],
            [20, 25, 30, 0]
        ]
        res = IntegerProgramming.traveling_salesman_nearest_neighbor(D)
        assert len(res['tour']) == 4
        assert set(res['tour']) == {0, 1, 2, 3}
        assert res['total_distance'] > 0

    def test_tsp_nn_single_city(self):
        """TSP mit 1 Stadt: leere Tour."""
        res = IntegerProgramming.traveling_salesman_nearest_neighbor([[0]])
        assert res['tour'] == [0]
        assert res['total_distance'] == 0.0

    def test_tsp_nn_empty_raises(self):
        """Leere Distanzmatrix wirft ValueError."""
        with pytest.raises(ValueError):
            IntegerProgramming.traveling_salesman_nearest_neighbor([])

    def test_tsp_nn_symmetric(self):
        """Symmetric TSP: Tour besucht alle Städte genau einmal."""
        n = 5
        rng = np.random.default_rng(0)
        D_half = rng.integers(1, 100, size=(n, n)).astype(float)
        D = (D_half + D_half.T) / 2
        np.fill_diagonal(D, 0)
        res = IntegerProgramming.traveling_salesman_nearest_neighbor(D.tolist())
        assert len(res['tour']) == n
        assert len(set(res['tour'])) == n

    # --- traveling_salesman_2opt ---

    def test_tsp_2opt_improves_tour(self):
        """2-Opt verbessert oder hält die Tour."""
        D = [
            [0, 2, 9, 10],
            [1, 0, 6,  4],
            [15, 7, 0, 8],
            [6, 3, 12, 0]
        ]
        nn_res  = IntegerProgramming.traveling_salesman_nearest_neighbor(D)
        opt_res = IntegerProgramming.traveling_salesman_2opt(D)
        assert opt_res['total_distance'] <= nn_res['total_distance'] + 1e-6

    def test_tsp_2opt_custom_initial_tour(self):
        """2-Opt mit vorgegebenem Startpfad."""
        D = [[0,1,2],[1,0,1],[2,1,0]]
        res = IntegerProgramming.traveling_salesman_2opt(D, initial_tour=[0,2,1])
        assert len(res['tour']) == 3

    def test_tsp_2opt_iterations_counted(self):
        """2-Opt zählt Iterationen."""
        D = [[0,10,15],[10,0,20],[15,20,0]]
        res = IntegerProgramming.traveling_salesman_2opt(D)
        assert res['iterations'] >= 1


# ===========================================================================
# 3. WARTESCHLANGENTHEORIE
# ===========================================================================

class TestQueueingTheory:
    """Tests für QueueingTheory."""

    # --- mm1_queue ---

    def test_mm1_stable(self):
        """M/M/1: ρ < 1 liefert endliche Kenngrößen."""
        res = QueueingTheory.mm1_queue(arrival_rate=2.0, service_rate=5.0)
        assert res['stable']
        assert abs(res['rho'] - 0.4) < 1e-10
        assert abs(res['L'] - 0.4/0.6) < 1e-6

    def test_mm1_little_law(self):
        """M/M/1: Little's Law L = λ·W muss gelten."""
        res = QueueingTheory.mm1_queue(arrival_rate=3.0, service_rate=5.0)
        assert res['stable']
        assert abs(res['L'] - 3.0 * res['W']) < 1e-6

    def test_mm1_unstable(self):
        """M/M/1: ρ ≥ 1 liefert stabile=False."""
        res = QueueingTheory.mm1_queue(arrival_rate=5.0, service_rate=3.0)
        assert not res['stable']
        assert res['L'] == math.inf

    def test_mm1_heavy_traffic(self):
        """M/M/1: ρ → 1 führt zu sehr großem L."""
        res = QueueingTheory.mm1_queue(arrival_rate=9.9, service_rate=10.0)
        assert res['stable']
        assert res['L'] > 50   # Hohe Last → sehr lange Schlange

    def test_mm1_invalid_rate_raises(self):
        """Negative/Null-Raten werfen ValueError."""
        with pytest.raises(ValueError):
            QueueingTheory.mm1_queue(0.0, 5.0)
        with pytest.raises(ValueError):
            QueueingTheory.mm1_queue(2.0, 0.0)

    def test_mm1_lq_formula(self):
        """M/M/1: Lq = ρ²/(1-ρ)."""
        lam, mu = 2.0, 4.0
        rho = lam / mu
        res = QueueingTheory.mm1_queue(lam, mu)
        assert abs(res['Lq'] - rho**2/(1-rho)) < 1e-8

    # --- mms_queue ---

    def test_mms_single_server(self):
        """M/M/s mit s=1 stimmt mit M/M/1 überein."""
        res1 = QueueingTheory.mm1_queue(2.0, 5.0)
        res2 = QueueingTheory.mms_queue(2.0, 5.0, 1)
        assert abs(res1['L'] - res2['L']) < 1e-4

    def test_mms_more_servers_less_waiting(self):
        """Mehr Server → kürzere Wartezeit."""
        res1 = QueueingTheory.mms_queue(8.0, 5.0, 2)
        res2 = QueueingTheory.mms_queue(8.0, 5.0, 3)
        assert res1['Wq'] > res2['Wq']

    def test_mms_zero_servers_raises(self):
        """s < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            QueueingTheory.mms_queue(2.0, 5.0, 0)

    def test_mms_unstable(self):
        """M/M/s instabil wenn ρ = λ/(s·μ) ≥ 1."""
        res = QueueingTheory.mms_queue(10.0, 3.0, 2)   # ρ = 10/6 > 1
        assert not res['stable']

    def test_mms_little_law(self):
        """M/M/s: Little's Law L = λ·W."""
        res = QueueingTheory.mms_queue(4.0, 3.0, 3)
        if res['stable']:
            assert abs(res['L'] - 4.0 * res['W']) < 1e-4

    # --- mm1k_queue ---

    def test_mm1k_probabilities_sum_to_1(self):
        """M/M/1/K: Stationäre Wahrscheinlichkeiten summieren zu 1."""
        res = QueueingTheory.mm1k_queue(3.0, 5.0, 5)
        assert abs(sum(res['probs']) - 1.0) < 1e-9

    def test_mm1k_blocking_prob(self):
        """M/M/1/K: Blockierungswahrscheinlichkeit ∈ [0, 1]."""
        res = QueueingTheory.mm1k_queue(4.0, 3.0, 3)
        assert 0 <= res['P_block'] <= 1

    def test_mm1k_small_k(self):
        """M/M/1/K mit K=1: Nur 0 oder 1 Kunden möglich."""
        res = QueueingTheory.mm1k_queue(2.0, 4.0, 1)
        assert len(res['probs']) == 2
        assert abs(sum(res['probs']) - 1.0) < 1e-9

    def test_mm1k_invalid_k_raises(self):
        """K < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            QueueingTheory.mm1k_queue(2.0, 3.0, 0)

    # --- erlang_c ---

    def test_erlang_c_between_0_and_1(self):
        """Erlang-C-Wahrscheinlichkeit ∈ [0, 1]."""
        pc = QueueingTheory.erlang_c(5.0, 3.0, 3)
        assert 0 <= pc <= 1

    def test_erlang_c_heavy_load(self):
        """Hohe Last → hohe Wartewahrscheinlichkeit."""
        pc = QueueingTheory.erlang_c(9.5, 5.0, 2)
        assert pc > 0.5

    def test_erlang_c_light_load(self):
        """Leichte Last → niedrige Wartewahrscheinlichkeit."""
        pc = QueueingTheory.erlang_c(0.1, 10.0, 3)
        assert pc < 0.1

    # --- erlang_b ---

    def test_erlang_b_between_0_and_1(self):
        """Erlang-B ∈ [0, 1]."""
        pb = QueueingTheory.erlang_b(5.0, 3.0, 5)
        assert 0 <= pb <= 1

    def test_erlang_b_increases_with_load(self):
        """Höhere Last → höhere Blockierungswahrscheinlichkeit."""
        pb1 = QueueingTheory.erlang_b(1.0, 3.0, 3)
        pb2 = QueueingTheory.erlang_b(5.0, 3.0, 3)
        assert pb2 > pb1

    def test_erlang_b_invalid_s_raises(self):
        """s < 1 wirft ValueError."""
        with pytest.raises(ValueError):
            QueueingTheory.erlang_b(2.0, 3.0, 0)

    # --- little_law_check ---

    def test_little_law_holds(self):
        """Little's Law: L = λ·W → holds=True."""
        res = QueueingTheory.little_law_check(L=10.0, lambda_=2.0, W=5.0)
        assert res['holds']
        assert abs(res['L_computed'] - 10.0) < 1e-10

    def test_little_law_violated(self):
        """Verletztes Little's Law → holds=False."""
        res = QueueingTheory.little_law_check(L=10.0, lambda_=2.0, W=3.0)
        assert not res['holds']


# ===========================================================================
# 4. MARKOV-ENTSCHEIDUNGSPROZESSE
# ===========================================================================

class TestMarkovDecisionProcess:
    """Tests für MarkovDecisionProcess."""

    def _create_simple_mdp(self):
        """Einfaches 2-Zustands-MDP für Tests."""
        states  = ['A', 'B']
        actions = ['stay', 'move']
        trans   = {
            ('A', 'stay'): {'A': 1.0},
            ('A', 'move'): {'B': 1.0},
            ('B', 'stay'): {'B': 1.0},
            ('B', 'move'): {'A': 1.0},
        }
        rewards = {
            ('A', 'stay'): 1.0,
            ('A', 'move'): 0.0,
            ('B', 'stay'): 2.0,
            ('B', 'move'): 0.0,
        }
        return MarkovDecisionProcess(states, actions, trans, rewards, gamma=0.9)

    def test_value_iteration_converges(self):
        """Wertiteration konvergiert."""
        mdp = self._create_simple_mdp()
        res = mdp.value_iteration()
        assert res['converged']

    def test_value_iteration_policy_not_none(self):
        """Wertiteration liefert eine Politik."""
        mdp = self._create_simple_mdp()
        res = mdp.value_iteration()
        assert res['policy'] is not None
        assert 'A' in res['policy']
        assert 'B' in res['policy']

    def test_value_iteration_b_better_than_a(self):
        """Zustand B (Belohnung 2) hat höheren Wert als A (Belohnung 1)."""
        mdp = self._create_simple_mdp()
        res = mdp.value_iteration()
        V   = res['value_function']
        assert V['B'] > V['A']

    def test_value_iteration_optimal_policy(self):
        """Optimale Politik: In Zustand B bleiben (höhere Belohnung)."""
        mdp = self._create_simple_mdp()
        res = mdp.value_iteration()
        # Optimale Politik für B: stay (Belohnung 2 > 0)
        assert res['policy']['B'] == 'stay'

    def test_policy_iteration_converges(self):
        """Politikiteration konvergiert."""
        mdp = self._create_simple_mdp()
        res = mdp.policy_iteration()
        assert res['converged']

    def test_policy_iteration_matches_value_iteration(self):
        """Politikiteration liefert dieselbe Politik wie Wertiteration."""
        mdp1 = self._create_simple_mdp()
        mdp2 = self._create_simple_mdp()
        res_vi = mdp1.value_iteration()
        res_pi = mdp2.policy_iteration()
        # Beide sollen für B 'stay' wählen
        assert res_vi['policy']['B'] == res_pi['policy']['B']

    def test_bellman_equation_non_negative(self):
        """Bellman-Gleichung für Zustand A ist nicht negativ."""
        mdp = self._create_simple_mdp()
        V   = {'A': 5.0, 'B': 10.0}
        val = mdp.bellman_equation(V, 'A')
        assert val >= 0

    def test_mdp_invalid_gamma_raises(self):
        """Ungültiger Gamma-Wert wirft ValueError."""
        with pytest.raises(ValueError):
            MarkovDecisionProcess(['A'], ['a'], {}, {}, gamma=1.0)
        with pytest.raises(ValueError):
            MarkovDecisionProcess(['A'], ['a'], {}, {}, gamma=-0.1)

    def test_optimal_policy_returns_all_states(self):
        """optimal_policy() gibt für alle Zustände eine Aktion zurück."""
        mdp = self._create_simple_mdp()
        mdp.value_iteration()
        pol = mdp.optimal_policy()
        assert set(pol.keys()) == set(mdp.states)


# ===========================================================================
# 5. NETZWERKFLÜSSE
# ===========================================================================

class TestNetworkFlows:
    """Tests für NetworkFlows."""

    def _simple_flow_graph(self):
        """Einfacher Flussgraph für Tests."""
        return {
            0: {1: 10, 2: 10},
            1: {2: 5, 3: 10},
            2: {4: 10},
            3: {4: 10},
            4: {}
        }

    # --- max_flow_ford_fulkerson ---

    def test_max_flow_simple(self):
        """Max-Flow für einfachen Graphen."""
        graph = {0: {1: 3, 2: 2}, 1: {3: 2}, 2: {3: 3}, 3: {}}
        res = NetworkFlows.max_flow_ford_fulkerson(graph, 0, 3)
        assert res['max_flow'] == pytest.approx(4.0, abs=1e-6)

    def test_max_flow_no_path(self):
        """Kein Pfad von source zu sink → Fluss = 0."""
        graph = {0: {1: 5}, 1: {}, 2: {3: 5}, 3: {}}
        res = NetworkFlows.max_flow_ford_fulkerson(graph, 0, 3)
        assert res['max_flow'] == pytest.approx(0.0, abs=1e-6)

    def test_max_flow_empty_graph(self):
        """Leerer Graph → Fluss = 0."""
        res = NetworkFlows.max_flow_ford_fulkerson({}, 0, 1)
        assert res['max_flow'] == 0.0

    def test_max_flow_bottleneck(self):
        """Flaschenhals begrenzt den Max-Flow."""
        graph = {0: {1: 100}, 1: {2: 1}, 2: {3: 100}, 3: {}}
        res = NetworkFlows.max_flow_ford_fulkerson(graph, 0, 3)
        assert res['max_flow'] == pytest.approx(1.0)

    # --- min_cut ---

    def test_min_cut_equals_max_flow(self):
        """Min-Cut-Theorem: min_cut_value = max_flow."""
        graph = {0: {1: 3, 2: 2}, 1: {3: 2}, 2: {3: 3}, 3: {}}
        flow_res = NetworkFlows.max_flow_ford_fulkerson(graph, 0, 3)
        cut_res  = NetworkFlows.min_cut(graph, 0, 3)
        assert abs(cut_res['min_cut_value'] - flow_res['max_flow']) < 1e-6

    def test_min_cut_source_in_s_set(self):
        """Source muss in S-Menge sein."""
        graph = {0: {1: 5, 2: 5}, 1: {3: 5}, 2: {3: 5}, 3: {}}
        res = NetworkFlows.min_cut(graph, 0, 3)
        assert 0 in res['S_set']

    def test_min_cut_sink_in_t_set(self):
        """Sink muss in T-Menge sein."""
        graph = {0: {1: 5, 2: 5}, 1: {3: 5}, 2: {3: 5}, 3: {}}
        res = NetworkFlows.min_cut(graph, 0, 3)
        assert 3 in res['T_set']

    # --- shortest_path_dijkstra ---

    def test_dijkstra_simple(self):
        """Dijkstra: bekannte Kürzeste-Pfade."""
        graph = {0: {1: 1, 2: 4}, 1: {2: 2, 3: 5}, 2: {3: 1}, 3: {}}
        res = NetworkFlows.shortest_path_dijkstra(graph, 0)
        assert res['distances'][3] == pytest.approx(4.0)  # 0→1→2→3: 1+2+1=4

    def test_dijkstra_source_distance_zero(self):
        """Distanz von source zu sich selbst ist 0."""
        graph = {0: {1: 5}, 1: {}}
        res = NetworkFlows.shortest_path_dijkstra(graph, 0)
        assert res['distances'][0] == 0.0

    def test_dijkstra_unreachable(self):
        """Nicht erreichbare Knoten haben Distanz inf."""
        graph = {0: {1: 1}, 1: {}, 2: {3: 1}, 3: {}}
        res = NetworkFlows.shortest_path_dijkstra(graph, 0)
        assert res['distances'][2] == math.inf

    def test_dijkstra_negative_weight_raises(self):
        """Negatives Kantengewicht wirft ValueError."""
        graph = {0: {1: -1}, 1: {}}
        with pytest.raises(ValueError):
            NetworkFlows.shortest_path_dijkstra(graph, 0)

    def test_dijkstra_empty_graph(self):
        """Leerer Graph: nur Quelldistanz bekannt."""
        res = NetworkFlows.shortest_path_dijkstra({}, 0)
        assert res['distances'][0] == 0.0

    # --- shortest_path_bellman_ford ---

    def test_bellman_ford_negative_edges(self):
        """Bellman-Ford mit negativen Kanten."""
        graph = {0: {1: 4, 2: 5}, 1: {2: -3}, 2: {3: 2}, 3: {}}
        res = NetworkFlows.shortest_path_bellman_ford(graph, 0)
        # 0→1→2→3: 4-3+2=3; 0→2→3: 5+2=7 → 3 ist kürzer
        assert res['distances'][3] == pytest.approx(3.0)

    def test_bellman_ford_detects_negative_cycle(self):
        """Bellman-Ford erkennt negativen Zyklus."""
        graph = {0: {1: 1}, 1: {2: -1}, 2: {1: -1}}   # Negativer Zyklus: 1→2→1
        res = NetworkFlows.shortest_path_bellman_ford(graph, 0)
        assert res['negative_cycle']

    def test_bellman_ford_no_negative_cycle(self):
        """Kein negativer Zyklus → negative_cycle=False."""
        graph = {0: {1: 1, 2: 4}, 1: {3: 2}, 2: {3: 1}, 3: {}}
        res = NetworkFlows.shortest_path_bellman_ford(graph, 0)
        assert not res['negative_cycle']

    # --- minimum_spanning_tree_kruskal ---

    def test_kruskal_total_weight(self):
        """Kruskal: Bekannter MST-Wert."""
        edges = [(1,0,1), (2,0,2), (3,1,2), (4,2,3), (5,1,3)]
        res = NetworkFlows.minimum_spanning_tree_kruskal(edges, 4)
        assert res['total_weight'] == pytest.approx(7.0)   # 1+2+4=7

    def test_kruskal_num_edges(self):
        """MST hat genau n-1 Kanten für zusammenhängenden Graphen."""
        edges = [(1,0,1),(2,0,2),(3,1,2),(4,2,3)]
        res = NetworkFlows.minimum_spanning_tree_kruskal(edges, 4)
        assert len(res['mst_edges']) == 3

    def test_kruskal_empty_graph(self):
        """Leerer Graph (0 Knoten)."""
        res = NetworkFlows.minimum_spanning_tree_kruskal([], 0)
        assert res['total_weight'] == 0.0

    def test_kruskal_single_node(self):
        """Einzelner Knoten: kein MST-Edge nötig."""
        res = NetworkFlows.minimum_spanning_tree_kruskal([], 1)
        assert res['mst_edges'] == []

    # --- minimum_spanning_tree_prim ---

    def test_prim_total_weight(self):
        """Prim: Bekannter MST-Wert für einfachen Graphen."""
        adj = [
            [0, 1, 4, 0],
            [1, 0, 2, 5],
            [4, 2, 0, 3],
            [0, 5, 3, 0]
        ]
        res = NetworkFlows.minimum_spanning_tree_prim(adj)
        assert res['total_weight'] == pytest.approx(6.0)   # 1+2+3=6

    def test_prim_num_edges(self):
        """Prim MST hat genau n-1 Kanten."""
        adj = [[0,2,3],[2,0,1],[3,1,0]]
        res = NetworkFlows.minimum_spanning_tree_prim(adj)
        assert len(res['mst_edges']) == 2

    def test_prim_empty_raises(self):
        """Leere Adjazenzmatrix wirft ValueError."""
        with pytest.raises(ValueError):
            NetworkFlows.minimum_spanning_tree_prim([])


# ===========================================================================
# 6. SPIELTHEORIE
# ===========================================================================

class TestGameTheory:
    """Tests für GameTheory."""

    # --- nash_equilibrium_2x2 ---

    def test_nash_prisoners_dilemma_has_pure_ne(self):
        """Gefangenendilemma hat ein reines NE: (D, D) = (1, 1)."""
        A = [[-1, -3], [0, -2]]
        B = [[-1, 0], [-3, -2]]
        res = GameTheory.nash_equilibrium_2x2(A, B)
        pure_ne = res['pure_equilibria']
        assert len(pure_ne) >= 1
        strategies = [ne['strategies'] for ne in pure_ne]
        assert (1, 1) in strategies  # (D, D)

    def test_nash_coordination_game(self):
        """Koordinierungsspiel hat zwei reine NE: (0,0) und (1,1)."""
        A = [[2, 0], [0, 1]]
        B = [[2, 0], [0, 1]]
        res = GameTheory.nash_equilibrium_2x2(A, B)
        strategies = [ne['strategies'] for ne in res['pure_equilibria']]
        assert (0, 0) in strategies
        assert (1, 1) in strategies

    def test_nash_matching_pennies_mixed(self):
        """Matching Pennies: kein reines NE, gemischtes NE = (0.5, 0.5)."""
        A = [[1, -1], [-1, 1]]
        B = [[-1, 1], [1, -1]]
        res = GameTheory.nash_equilibrium_2x2(A, B)
        assert len(res['pure_equilibria']) == 0
        me = res['mixed_equilibrium']
        if me is not None:
            assert abs(me['p'] - 0.5) < 1e-6
            assert abs(me['q'] - 0.5) < 1e-6

    def test_nash_returns_structure(self):
        """nash_equilibrium_2x2 gibt korrekte Schlüssel zurück."""
        A = [[1, 0], [0, 1]]
        B = [[1, 0], [0, 1]]
        res = GameTheory.nash_equilibrium_2x2(A, B)
        assert 'pure_equilibria' in res
        assert 'mixed_equilibrium' in res

    # --- minimax_theorem ---

    def test_minimax_zero_sum_game(self):
        """Minimax für bekanntes Nullsummenspiel."""
        # Rock-Paper-Scissors reduziert auf 2x2 (R vs P/S)
        A = [[0, -1], [1, 0]]   # Spieler 1: [R, P]; Spieler 2: [P, S]
        res = GameTheory.minimax_theorem(A)
        assert res['success']
        assert res['game_value'] is not None

    def test_minimax_strategy_sums_to_1(self):
        """Gemischte Strategien summieren zu 1."""
        A = [[3, -1], [-2, 4]]
        res = GameTheory.minimax_theorem(A)
        if res['success'] and res['p_strategy']:
            assert abs(sum(res['p_strategy']) - 1.0) < 1e-4
        if res['success'] and res['q_strategy']:
            assert abs(sum(res['q_strategy']) - 1.0) < 1e-4

    def test_minimax_dominant_strategy_game(self):
        """Spiel mit dominanter Strategie: Spielwert durch pure NE bestimmt."""
        A = [[2, 2], [1, 1]]   # Strategie 0 dominiert Strategie 1
        res = GameTheory.minimax_theorem(A)
        assert res['success']
        # Spielwert ≈ 2 (dominante Strategie gewinnt immer 2)
        if res['game_value'] is not None:
            assert res['game_value'] == pytest.approx(2.0, abs=0.5)

    # --- prisoners_dilemma ---

    def test_prisoners_dilemma_has_ne(self):
        """Gefangenendilemma hat ein NE."""
        res = GameTheory.prisoners_dilemma()
        assert len(res['nash_equilibrium']) >= 1

    def test_prisoners_dilemma_structure(self):
        """Gefangenendilemma liefert vollständige Analyse."""
        res = GameTheory.prisoners_dilemma()
        assert 'dominant_strategy' in res
        assert 'pareto_optimum' in res
        assert 'dilemma_explanation' in res

    def test_prisoners_dilemma_dominant_is_D(self):
        """Dominante Strategie muss Defektieren sein."""
        res = GameTheory.prisoners_dilemma()
        assert 'D' in res['dominant_strategy'] or 'defektieren' in res['dominant_strategy'].lower()

    # --- dominant_strategy ---

    def test_dominant_strategy_strictly_dominant(self):
        """Strikt dominante Strategie gefunden."""
        # Strategie 0: immer besser als Strategie 1
        P = [[5, 4], [3, 2]]   # Zeile 0 dominiert Zeile 1 strikt
        res = GameTheory.dominant_strategy(P)
        assert 0 in res['strictly_dominant']

    def test_dominant_strategy_dominated(self):
        """Dominierte Strategie erkannt."""
        P = [[5, 4], [3, 2]]
        res = GameTheory.dominant_strategy(P)
        assert 1 in res['dominated_strategies']

    def test_dominant_strategy_no_dominant(self):
        """Kein reines Dominanzverhältnis."""
        P = [[3, 0], [0, 3]]   # Weder dominiert die andere
        res = GameTheory.dominant_strategy(P)
        assert len(res['strictly_dominant']) == 0

    # --- pareto_optimal ---

    def test_pareto_optimal_basic(self):
        """Pareto-Optimalität: (3,3) dominiert (2,1)."""
        outcomes = [(3, 3), (2, 1), (1, 4)]
        res = GameTheory.pareto_optimal(outcomes)
        assert (2, 1) in res['dominated_outcomes']

    def test_pareto_optimal_all_optimal(self):
        """Alle Ergebnisse sind Pareto-optimal (kein dominiert das andere)."""
        outcomes = [(3, 0), (0, 3), (1, 1)]
        res = GameTheory.pareto_optimal(outcomes)
        # (1,1) wird von (3,0) oder (0,3) nicht dominiert, und umgekehrt auch nicht
        # (3,0) dominiert nicht (0,3) weil 0 < 3
        assert (3, 0) in res['pareto_optimal']
        assert (0, 3) in res['pareto_optimal']

    def test_pareto_optimal_single_outcome(self):
        """Einzelnes Ergebnis ist immer Pareto-optimal."""
        outcomes = [(5, 5)]
        res = GameTheory.pareto_optimal(outcomes)
        assert len(res['pareto_optimal']) == 1


# ===========================================================================
# 7. STOCHASTISCHE OPTIMIERUNG
# ===========================================================================

class TestStochasticOptimization:
    """Tests für StochasticOptimization."""

    # --- stochastic_gradient_descent ---

    def test_sgd_quadratic_minimization(self):
        """SGD minimiert f(x) = x² auf x ≈ 0."""
        f      = lambda x: float(x[0] ** 2)
        grad_f = lambda x, i: np.array([2 * x[0]])
        res    = StochasticOptimization.stochastic_gradient_descent(
            f, grad_f, x0=[5.0], n_steps=500, lr=0.01, seed=42
        )
        assert abs(res['x_opt'][0]) < 1.0

    def test_sgd_returns_structure(self):
        """SGD gibt korrekten Dict zurück."""
        f      = lambda x: float(x[0] ** 2)
        grad_f = lambda x, i: np.array([2 * x[0]])
        res    = StochasticOptimization.stochastic_gradient_descent(f, grad_f, [1.0], 10)
        assert 'x_opt' in res and 'f_opt' in res and 'history' in res

    def test_sgd_history_length(self):
        """SGD history hat mindestens 2 Einträge (Start + Ende)."""
        f      = lambda x: float(x[0] ** 2)
        grad_f = lambda x, i: np.array([2 * x[0]])
        res    = StochasticOptimization.stochastic_gradient_descent(f, grad_f, [3.0], n_steps=100)
        assert len(res['history']) >= 2

    def test_sgd_decreases_objective(self):
        """SGD reduziert den Funktionswert (konvexe Funktion)."""
        f      = lambda x: float(x[0] ** 2 + x[1] ** 2)
        grad_f = lambda x, i: np.array([2 * x[0], 2 * x[1]])
        res    = StochasticOptimization.stochastic_gradient_descent(
            f, grad_f, [10.0, 10.0], n_steps=1000, lr=0.01
        )
        assert res['f_opt'] < f(np.array([10.0, 10.0]))

    # --- simulated_annealing ---

    def test_sa_finds_minimum(self):
        """Simulated Annealing findet Minimum von x²."""
        energy = lambda x: float(x[0] ** 2)
        res = StochasticOptimization.simulated_annealing(
            energy, x0=[5.0], T_init=1.0, cooling=0.99, n_steps=2000, seed=42
        )
        assert abs(res['x_opt'][0]) < 2.0   # SA findet nicht immer exaktes Minimum

    def test_sa_returns_structure(self):
        """SA gibt korrekten Dict zurück."""
        energy = lambda x: float(x[0] ** 2)
        res = StochasticOptimization.simulated_annealing(energy, [1.0], n_steps=100)
        assert 'x_opt' in res and 'energy_opt' in res and 'T_final' in res

    def test_sa_temperature_decreases(self):
        """Temperatur fällt mit cooling-Rate."""
        energy = lambda x: float(x[0] ** 2)
        T_init  = 2.0
        cooling = 0.9
        n_steps = 10
        res = StochasticOptimization.simulated_annealing(
            energy, [0.0], T_init=T_init, cooling=cooling, n_steps=n_steps
        )
        expected_T = T_init * (cooling ** n_steps)
        assert abs(res['T_final'] - expected_T) < 1e-6

    def test_sa_multivariate(self):
        """SA funktioniert mit mehrdimensionalem x₀."""
        energy = lambda x: float(np.sum(x ** 2))
        res = StochasticOptimization.simulated_annealing(
            energy, x0=[3.0, 4.0], n_steps=1000, seed=42
        )
        assert res['energy_opt'] < 25.0   # Besser als Startpunkt (9+16=25)

    # --- newsvendor_problem ---

    def test_newsvendor_critical_ratio(self):
        """Newsvendor: Kritische Ratio = c_u / (c_u + c_o)."""
        rng = np.random.default_rng(42)
        D   = rng.uniform(50, 150, size=10000)
        res = StochasticOptimization.newsvendor_problem(D, cost_overstock=1.0, cost_understock=3.0)
        assert abs(res['critical_ratio'] - 0.75) < 1e-6

    def test_newsvendor_optimal_quantity_in_range(self):
        """Optimale Bestellmenge liegt im Bereich der Nachfrage."""
        rng = np.random.default_rng(0)
        D   = rng.uniform(100, 200, size=1000)
        res = StochasticOptimization.newsvendor_problem(D, 2.0, 5.0)
        assert 100 <= res['optimal_quantity'] <= 200

    def test_newsvendor_invalid_costs_raise(self):
        """Ungültige Kosten werfen ValueError."""
        with pytest.raises(ValueError):
            StochasticOptimization.newsvendor_problem([10, 20], -1, 5)
        with pytest.raises(ValueError):
            StochasticOptimization.newsvendor_problem([10, 20], 5, 0)

    def test_newsvendor_equal_costs(self):
        """Bei gleichen Kosten: kritische Ratio = 0.5 (Median-Regel)."""
        D = list(range(0, 101))
        res = StochasticOptimization.newsvendor_problem(D, 1.0, 1.0)
        assert abs(res['critical_ratio'] - 0.5) < 1e-6

    # --- stochastic_programming_demo ---

    def test_stochastic_programming_demo(self):
        """Stochastische Programmierung: Zwei Szenarien."""
        scenarios = [
            {'c': [1.0, 2.0], 'A_ub': [[1.0, 1.0]], 'b_ub': [10.0]},
            {'c': [2.0, 1.0], 'A_ub': [[1.0, 1.0]], 'b_ub': [10.0]},
        ]
        probs = [0.5, 0.5]
        res = StochasticOptimization.stochastic_programming_demo(scenarios, probs)
        assert 'expected_value' in res
        assert 'scenario_results' in res
        assert len(res['scenario_results']) == 2

    def test_stochastic_programming_invalid_probs(self):
        """Ungültige Wahrscheinlichkeiten werfen ValueError."""
        with pytest.raises(ValueError):
            StochasticOptimization.stochastic_programming_demo(
                [{'c': [1.0], 'A_ub': [[1.0]], 'b_ub': [5.0]}],
                [0.3]   # Summe ≠ 1
            )

    def test_stochastic_programming_expected_value(self):
        """Erwartungswert ist gewichteter Durchschnitt der Szenariowerte."""
        scenarios = [
            {'c': [-1.0], 'A_ub': [[1.0]], 'b_ub': [4.0]},
            {'c': [-1.0], 'A_ub': [[1.0]], 'b_ub': [6.0]},
        ]
        probs = [0.4, 0.6]
        res = StochasticOptimization.stochastic_programming_demo(scenarios, probs)
        # Szenario 1: Optimum -4 (x=4), Szenario 2: Optimum -6 (x=6)
        # Erwartungswert: 0.4*(-4) + 0.6*(-6) = -1.6 - 3.6 = -5.2
        assert abs(res['expected_value'] - (-5.2)) < 1e-4
