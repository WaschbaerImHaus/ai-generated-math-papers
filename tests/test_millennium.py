"""
@file test_millennium.py
@brief Tests für das millennium_problems-Modul.
@description
    Umfassende Tests für alle implementierten Werkzeuge zu den
    7 Millennium-Problemen des Clay Mathematics Institute.
    Testet Korrektheit, Konsistenz und Edge-Cases.
@author Kurt Ingwer
@lastModified 2026-03-10
"""

import pytest
import math
import numpy as np
import sys
import os

# Projektpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from millennium_problems import (
    riemann_zeros_mpmath,
    verify_rh_critical_line,
    gram_law_check,
    zeta_zero_gaps,
    goldbach_verification_range,
    goldbach_weak_check,
    hardy_littlewood_goldbach_density,
    is_np_complete_3sat,
    traveling_salesman_brute,
    traveling_salesman_greedy,
    complexity_comparison,
    navier_stokes_2d_simple,
    check_navier_stokes_energy,
    yang_mills_action,
    lattice_gauge_theory_demo,
    elliptic_curve_points_mod_p,
    bsd_l_function_rank_estimate,
    millennium_problems_status,
)


# =============================================================================
# TESTS: RIEMANN-HYPOTHESE
# =============================================================================

class TestRiemannHypothesis:
    """Tests für die Riemann-Hypothese Werkzeuge."""

    def test_riemann_zeros_mpmath_count(self):
        """riemann_zeros_mpmath(5) gibt genau 5 Einträge zurück."""
        zeros = riemann_zeros_mpmath(5)
        assert len(zeros) == 5

    def test_riemann_zeros_mpmath_structure(self):
        """Jeder Eintrag hat die erwarteten Schlüssel."""
        zeros = riemann_zeros_mpmath(3)
        for z in zeros:
            assert 'n' in z
            assert 'zero' in z
            assert 're' in z
            assert 'im' in z
            assert 'on_line' in z

    def test_riemann_zeros_on_critical_line(self):
        """Alle berechneten Nullstellen liegen auf Re(s) = 1/2."""
        zeros = riemann_zeros_mpmath(5)
        for z in zeros:
            assert z['on_line'], f"Nullstelle {z['n']} nicht auf Re=1/2: re={z['re']}"
            assert abs(z['re'] - 0.5) < 1e-8

    def test_riemann_zeros_positive_imaginary(self):
        """Alle Imaginärteile der Nullstellen sind positiv."""
        zeros = riemann_zeros_mpmath(5)
        for z in zeros:
            assert z['im'] > 0

    def test_riemann_zeros_increasing(self):
        """Imaginärteile sind aufsteigend geordnet."""
        zeros = riemann_zeros_mpmath(5)
        im_parts = [z['im'] for z in zeros]
        assert im_parts == sorted(im_parts)

    def test_riemann_zeros_first_known(self):
        """Erste bekannte Nullstelle liegt bei t ≈ 14.134725..."""
        zeros = riemann_zeros_mpmath(1)
        assert abs(zeros[0]['im'] - 14.134725) < 0.001

    def test_riemann_zeros_indexing(self):
        """n-Felder sind korrekt nummeriert (1, 2, 3, ...)."""
        zeros = riemann_zeros_mpmath(5)
        for i, z in enumerate(zeros):
            assert z['n'] == i + 1

    def test_verify_rh_critical_line_basic(self):
        """verify_rh_critical_line(10) → all_on_line=True."""
        result = verify_rh_critical_line(10)
        assert result['all_on_line'] is True

    def test_verify_rh_critical_line_structure(self):
        """Rückgabe enthält erwartete Felder."""
        result = verify_rh_critical_line(5)
        assert 'verified' in result
        assert 'max_deviation' in result
        assert 'all_on_line' in result
        assert 'total' in result

    def test_verify_rh_100_zeros(self):
        """Erste 100 Nullstellen alle auf kritischer Geraden."""
        result = verify_rh_critical_line(100)
        assert result['all_on_line'] is True
        assert result['verified'] == 100

    def test_verify_rh_deviation_small(self):
        """Maximale Abweichung von 0.5 ist sehr klein."""
        result = verify_rh_critical_line(10)
        assert result['max_deviation'] < 1e-10

    def test_gram_law_check_structure(self):
        """gram_law_check(5) gibt dict mit korrekter Struktur zurück."""
        result = gram_law_check(5)
        assert isinstance(result, dict)
        assert 'gram_points' in result
        assert 'law_holds_count' in result
        assert 'total' in result
        assert 'success_rate' in result

    def test_gram_law_check_count(self):
        """gram_law_check(n) enthält genau n Gram-Punkte."""
        result = gram_law_check(5)
        assert result['total'] == 5
        assert len(result['gram_points']) == 5

    def test_gram_law_check_gram_points_positive(self):
        """Alle Gram-Punkte sind positive reelle Zahlen."""
        result = gram_law_check(5)
        for gp in result['gram_points']:
            assert gp['gram_point'] > 0

    def test_gram_law_check_success_rate_range(self):
        """Erfolgsquote liegt zwischen 0 und 1."""
        result = gram_law_check(5)
        assert 0.0 <= result['success_rate'] <= 1.0

    def test_zeta_zero_gaps_basic(self):
        """zeta_zero_gaps(10) gibt dict mit 'gaps' zurück."""
        result = zeta_zero_gaps(10)
        assert isinstance(result, dict)
        assert 'gaps' in result

    def test_zeta_zero_gaps_count(self):
        """Anzahl der Lücken ist n-1 bei n Nullstellen."""
        result = zeta_zero_gaps(10)
        assert len(result['gaps']) == 9  # 10 Nullstellen → 9 Lücken

    def test_zeta_zero_gaps_positive(self):
        """Alle Lücken zwischen Nullstellen sind positiv."""
        result = zeta_zero_gaps(10)
        for gap in result['gaps']:
            assert gap > 0

    def test_zeta_zero_gaps_statistics(self):
        """Statistische Kennwerte sind vorhanden und konsistent."""
        result = zeta_zero_gaps(20)
        assert result['mean'] > 0
        assert result['std'] >= 0
        assert result['min'] <= result['mean']
        assert result['max'] >= result['mean']


# =============================================================================
# TESTS: GOLDBACH-VERMUTUNG
# =============================================================================

class TestGoldbachConjecture:
    """Tests für die Goldbach-Vermutung Werkzeuge."""

    def test_goldbach_verification_range_small(self):
        """Goldbach-Vermutung für alle geraden n ≤ 100 verifiziert."""
        result = goldbach_verification_range(100)
        assert result['all_verified'] is True

    def test_goldbach_verification_range_medium(self):
        """Goldbach-Vermutung für alle geraden n ≤ 10000 verifiziert."""
        result = goldbach_verification_range(10000)
        assert result['all_verified'] is True
        assert result['verified_up_to'] == 10000

    def test_goldbach_verification_range_structure(self):
        """Rückgabe enthält erwartete Felder."""
        result = goldbach_verification_range(50)
        assert 'verified_up_to' in result
        assert 'all_verified' in result
        assert 'min_decompositions' in result

    def test_goldbach_weak_check_9(self):
        """Schwache Goldbach-Vermutung: 9 = p1 + p2 + p3."""
        result = goldbach_weak_check(9)
        assert result['verified'] is True
        assert result['decomposition'] is not None
        p1, p2, p3 = result['decomposition']
        assert p1 + p2 + p3 == 9

    def test_goldbach_weak_check_primes(self):
        """Alle drei Summanden sind Primzahlen."""
        from millennium_problems import _simple_sieve
        primes = set(_simple_sieve(100))
        for n in [7, 9, 11, 13, 15, 21, 27, 99]:
            result = goldbach_weak_check(n)
            assert result['verified'], f"Kein Beweis für n={n}"
            p1, p2, p3 = result['decomposition']
            assert p1 in primes, f"{p1} ist keine Primzahl"
            assert p2 in primes, f"{p2} ist keine Primzahl"
            assert p3 in primes, f"{p3} ist keine Primzahl"

    def test_goldbach_weak_check_invalid_even(self):
        """Schwache Vermutung gilt nur für ungerade Zahlen."""
        result = goldbach_weak_check(8)
        assert result['verified'] is False

    def test_goldbach_weak_check_too_small(self):
        """Für n ≤ 5 kein Beweis (Vermutung gilt erst für n > 5)."""
        result = goldbach_weak_check(5)
        assert result['verified'] is False

    def test_hardy_littlewood_density_structure(self):
        """hardy_littlewood_goldbach_density gibt korrekte Felder zurück."""
        result = hardy_littlewood_goldbach_density(100)
        assert 'estimate_r2' in result
        assert 'actual_r2' in result
        assert 'relative_error' in result
        assert 'twin_prime_constant_c2' in result

    def test_hardy_littlewood_density_positive(self):
        """Schätzung für n=100 ist positiv."""
        result = hardy_littlewood_goldbach_density(100)
        assert result['estimate_r2'] > 0
        assert result['actual_r2'] > 0

    def test_hardy_littlewood_c2_range(self):
        """Zwillingsprimzahl-Konstante liegt nahe bei 0.6601618."""
        result = hardy_littlewood_goldbach_density(100)
        c2 = result['twin_prime_constant_c2']
        assert 0.50 < c2 < 0.80

    def test_hardy_littlewood_invalid_input(self):
        """Ungerade Eingabe gibt Fehler zurück."""
        result = hardy_littlewood_goldbach_density(7)
        assert 'error' in result


# =============================================================================
# TESTS: P vs NP / KOMPLEXITÄTSTHEORIE
# =============================================================================

class TestPvsNP:
    """Tests für die P vs NP Komplexitätswerkzeuge."""

    def test_3sat_satisfiable_simple(self):
        """Einfache erfüllbare 3-SAT-Instanz."""
        result = is_np_complete_3sat([[1, 2, 3]], 3)
        assert result['satisfiable'] is True

    def test_3sat_assignment_valid(self):
        """Gefundene Belegung erfüllt tatsächlich alle Klauseln."""
        clauses = [[1, -2, 3], [-1, 2, 3]]
        result = is_np_complete_3sat(clauses, 3)
        assert result['satisfiable'] is True
        assignment = result['assignment']
        for clause in clauses:
            clause_val = any(
                assignment.get(abs(lit), False) == (lit > 0)
                for lit in clause
            )
            assert clause_val, f"Klausel {clause} nicht erfüllt"

    def test_3sat_unsatisfiable(self):
        """Unerfüllbare 3-SAT-Instanz (x1 ∧ ¬x1 Widerspruch)."""
        # (x1 ∨ x1 ∨ x1) ∧ (¬x1 ∨ ¬x1 ∨ ¬x1) → unerfüllbar
        clauses = [[1, 1, 1], [-1, -1, -1]]
        result = is_np_complete_3sat(clauses, 1)
        assert result['satisfiable'] is False
        assert result['assignment'] is None

    def test_3sat_steps_counted(self):
        """Schritte werden gezählt."""
        result = is_np_complete_3sat([[1, 2, 3]], 3)
        assert result['steps'] > 0

    def test_3sat_empty_clauses(self):
        """Leere Klauselliste ist trivial erfüllbar."""
        result = is_np_complete_3sat([], 3)
        assert result['satisfiable'] is True

    def test_tsp_brute_small(self):
        """TSP Brute-Force für 3 Städte."""
        dist = [[0, 1, 2], [1, 0, 1], [2, 1, 0]]
        result = traveling_salesman_brute(dist)
        assert 'tour' in result
        assert 'length' in result
        assert result['optimal'] is True
        # Tour enthält alle Städte (0, 1, 2)
        assert set(result['tour']) == {0, 1, 2}

    def test_tsp_brute_optimal_length(self):
        """Brute-Force findet optimale Tour."""
        dist = [[0, 10, 15, 20], [10, 0, 35, 25], [15, 35, 0, 30], [20, 25, 30, 0]]
        result = traveling_salesman_brute(dist)
        # Optimale Tour: 0→1→3→2→0 = 10+25+30+15 = 80
        assert result['length'] == 80.0

    def test_tsp_brute_too_large(self):
        """Brute-Force verweigert n > 12."""
        dist = [[0] * 13 for _ in range(13)]
        with pytest.raises(ValueError):
            traveling_salesman_brute(dist)

    def test_tsp_greedy_basic(self):
        """TSP Greedy für einfache Distanzmatrix."""
        dist = [[0, 1, 2], [1, 0, 1], [2, 1, 0]]
        result = traveling_salesman_greedy(dist)
        assert 'tour' in result
        assert 'length' in result
        assert result['optimal'] is False
        assert len(result['tour']) == 3
        assert set(result['tour']) == {0, 1, 2}

    def test_tsp_greedy_all_cities_visited(self):
        """Greedy-Tour besucht alle Städte genau einmal."""
        n = 5
        dist = [[abs(i - j) for j in range(n)] for i in range(n)]
        result = traveling_salesman_greedy(dist)
        assert len(result['tour']) == n
        assert len(set(result['tour'])) == n

    def test_tsp_greedy_length_positive(self):
        """Tourlänge ist positiv."""
        dist = [[0, 3, 5], [3, 0, 4], [5, 4, 0]]
        result = traveling_salesman_greedy(dist)
        assert result['length'] > 0

    def test_complexity_comparison_structure(self):
        """complexity_comparison gibt dict mit allen Komplexitäten zurück."""
        result = complexity_comparison(10)
        assert 'n' in result
        assert 'complexities' in result
        complexities = result['complexities']
        expected_keys = ['O(1)', 'O(log n)', 'O(n)', 'O(n log n)',
                         'O(n^2)', 'O(n^3)', 'O(2^n)', 'O(n!)']
        for key in expected_keys:
            assert key in complexities, f"Klasse {key} fehlt"

    def test_complexity_comparison_values(self):
        """Komplexitätswerte sind korrekt für n=10."""
        result = complexity_comparison(10)
        c = result['complexities']
        assert c['O(1)'] == 1
        assert c['O(n)'] == 10
        assert c['O(n^2)'] == 100
        assert abs(c['O(log n)'] - math.log2(10)) < 1e-10

    def test_complexity_comparison_ordering(self):
        """Komplexitäten sind aufsteigend geordnet."""
        result = complexity_comparison(10)
        c = result['complexities']
        # O(1) < O(log n) < O(n) < O(n²) < O(2^n) für n=10
        assert c['O(1)'] < c['O(log n)']
        assert c['O(log n)'] < c['O(n)']
        assert c['O(n)'] < c['O(n^2)']
        assert c['O(n^2)'] < c['O(2^n)']

    def test_complexity_comparison_n1(self):
        """Edge-Case: n=1."""
        result = complexity_comparison(1)
        assert result['n'] == 1
        assert result['complexities']['O(n)'] == 1


# =============================================================================
# TESTS: NAVIER-STOKES
# =============================================================================

class TestNavierStokes:
    """Tests für die Navier-Stokes Werkzeuge."""

    def test_navier_stokes_returns_dict(self):
        """navier_stokes_2d_simple gibt dict zurück."""
        result = navier_stokes_2d_simple(20, 20, T=0.1)
        assert isinstance(result, dict)

    def test_navier_stokes_fields_present(self):
        """Ergebnis enthält u, v, p, time, converged."""
        result = navier_stokes_2d_simple(20, 20, T=0.1)
        assert 'u' in result
        assert 'v' in result
        assert 'p' in result
        assert 'time' in result
        assert 'converged' in result

    def test_navier_stokes_grid_shape(self):
        """Ausgabe-Arrays haben korrekte Form."""
        nx, ny = 15, 20
        result = navier_stokes_2d_simple(nx, ny, T=0.05)
        assert result['u'].shape == (ny, nx)
        assert result['v'].shape == (ny, nx)
        assert result['p'].shape == (ny, nx)

    def test_navier_stokes_converged(self):
        """Simulation konvergiert für moderate Parameter."""
        result = navier_stokes_2d_simple(20, 20, dt=0.001, T=0.01, nu=0.1)
        assert result['converged'] is True

    def test_navier_stokes_boundary_conditions(self):
        """Randbedingungen: u=0 an Seitenwänden und Boden."""
        result = navier_stokes_2d_simple(20, 20, T=0.05)
        u = result['u']
        # Boden: u[0, :] = 0
        assert np.allclose(u[0, :], 0.0, atol=1e-12)
        # Links und Rechts: u[:, 0] = 0, u[:, -1] = 0
        assert np.allclose(u[:, 0], 0.0, atol=1e-12)
        assert np.allclose(u[:, -1], 0.0, atol=1e-12)

    def test_navier_stokes_reynolds_number(self):
        """Reynolds-Zahl wird berechnet."""
        result = navier_stokes_2d_simple(20, 20, T=0.05, nu=0.01)
        assert 'reynolds_number' in result
        assert abs(result['reynolds_number'] - 100.0) < 1e-6

    def test_check_energy_basic(self):
        """check_navier_stokes_energy gibt Energiewert zurück."""
        u = np.ones((10, 10)) * 0.5
        v = np.zeros((10, 10))
        result = check_navier_stokes_energy(u, v)
        assert 'kinetic_energy' in result
        assert result['kinetic_energy'] > 0

    def test_check_energy_zero_field(self):
        """Kinetische Energie ist 0 für Nullfeld."""
        u = np.zeros((10, 10))
        v = np.zeros((10, 10))
        result = check_navier_stokes_energy(u, v)
        assert abs(result['kinetic_energy']) < 1e-12

    def test_check_energy_bounded(self):
        """energy_bounded ist True für reguläre Felder."""
        u = np.random.rand(10, 10)
        v = np.random.rand(10, 10)
        result = check_navier_stokes_energy(u, v)
        assert result['energy_bounded'] is True

    def test_check_energy_structure(self):
        """Rückgabe enthält alle Felder."""
        u = np.ones((5, 5))
        v = np.ones((5, 5))
        result = check_navier_stokes_energy(u, v)
        assert 'kinetic_energy' in result
        assert 'max_velocity' in result
        assert 'energy_bounded' in result
        assert 'enstrophy' in result


# =============================================================================
# TESTS: YANG-MILLS
# =============================================================================

class TestYangMills:
    """Tests für die Yang-Mills Werkzeuge."""

    def test_yang_mills_action_positive(self):
        """Yang-Mills-Wirkung ist nicht-negativ."""
        F = np.array([[1.0, -0.5], [0.5, 0.0]])
        action = yang_mills_action(F)
        assert action >= 0

    def test_yang_mills_action_zero(self):
        """Nullfeld ergibt Wirkung 0."""
        F = np.zeros((4, 4))
        action = yang_mills_action(F)
        assert action == 0.0

    def test_yang_mills_action_coupling(self):
        """Wirkung skaliert mit 1/g²."""
        F = np.array([[1.0, 1.0], [1.0, 1.0]])
        s1 = yang_mills_action(F, coupling=1.0)
        s2 = yang_mills_action(F, coupling=2.0)
        # S ∝ 1/g², also S(g=2) = S(g=1)/4
        assert abs(s2 - s1 / 4.0) < 1e-12

    def test_yang_mills_action_zero_coupling(self):
        """g=0 wirft ValueError."""
        F = np.ones((3, 3))
        with pytest.raises(ValueError):
            yang_mills_action(F, coupling=0.0)

    def test_lattice_gauge_theory_structure(self):
        """lattice_gauge_theory_demo gibt dict mit wilson_loops zurück."""
        result = lattice_gauge_theory_demo(4, steps=10)
        assert 'wilson_loops' in result
        assert 'mass_gap_estimate' in result
        assert 'string_tension' in result

    def test_lattice_gauge_theory_wilson_loops(self):
        """Wilson-Loops liegen im erwarteten Bereich [-1, 1]."""
        result = lattice_gauge_theory_demo(4, steps=10)
        for wl in result['wilson_loops']:
            assert -1.1 <= wl['W'] <= 1.1

    def test_lattice_gauge_theory_mass_gap_non_negative(self):
        """Massenlücke-Schätzung ist nicht-negativ."""
        result = lattice_gauge_theory_demo(4, steps=10)
        assert result['mass_gap_estimate'] >= 0


# =============================================================================
# TESTS: BSD-VERMUTUNG
# =============================================================================

class TestBSDConjecture:
    """Tests für Birch-Swinnerton-Dyer Werkzeuge."""

    def test_elliptic_curve_points_mod_p_count(self):
        """Punkte auf y²=x³-x über F_5 werden gezählt."""
        # Kurve y²=x³-x (a=-1, b=0) über F_5
        points = elliptic_curve_points_mod_p(-1, 0, 5)
        # Enthält immer den Punkt im Unendlichen
        assert 'O' in points

    def test_elliptic_curve_points_hasse_bound(self):
        """Hasse-Weil-Schranke: |#E(F_p) - (p+1)| ≤ 2√p."""
        p = 7
        points = elliptic_curve_points_mod_p(1, 1, p)
        n = len(points)  # Inklusive Punkt im Unendlichen
        assert abs(n - (p + 1)) <= 2 * math.sqrt(p) + 1  # +1 für Rundungstoleranzen

    def test_elliptic_curve_singular_raises(self):
        """Singuläre Kurve (Diskriminante=0) wirft ValueError."""
        # y²=x³ ist singulär: a=0, b=0 → 4·0³+27·0²=0
        with pytest.raises(ValueError):
            elliptic_curve_points_mod_p(0, 0, 5)

    def test_bsd_rank_estimate_structure(self):
        """bsd_l_function_rank_estimate gibt korrekte Felder zurück."""
        from millennium_problems import _simple_sieve
        primes = _simple_sieve(50)
        result = bsd_l_function_rank_estimate(1, -1, primes)
        assert 'rank_estimate' in result
        assert 'log_product' in result
        assert 'curve' in result

    def test_bsd_rank_estimate_non_negative(self):
        """Rang-Schätzung ist nicht-negativ."""
        from millennium_problems import _simple_sieve
        primes = _simple_sieve(100)
        result = bsd_l_function_rank_estimate(-1, 0, primes)
        assert result['rank_estimate'] >= 0


# =============================================================================
# TESTS: STATUS-REPORT
# =============================================================================

class TestMillenniumStatus:
    """Tests für den Millennium-Probleme Status-Report."""

    def test_status_has_7_problems(self):
        """millennium_problems_status() hat genau 7 Einträge."""
        result = millennium_problems_status()
        assert len(result['problems']) == 7

    def test_status_keys_present(self):
        """Alle 7 Probleme sind vorhanden."""
        result = millennium_problems_status()
        expected_keys = [
            'P_vs_NP', 'Hodge_Conjecture', 'Poincare_Conjecture',
            'Riemann_Hypothesis', 'Yang_Mills', 'Navier_Stokes', 'BSD_Conjecture'
        ]
        for key in expected_keys:
            assert key in result['problems'], f"Problem {key} fehlt"

    def test_poincare_solved(self):
        """Poincaré-Vermutung ist als GELÖST markiert."""
        result = millennium_problems_status()
        poincare = result['problems']['Poincare_Conjecture']
        assert poincare['status'] == 'GELÖST'
        assert poincare['solved_by'] == 'Grigori Perelman'
        assert poincare['year_solved'] == 2003

    def test_other_problems_open(self):
        """Alle anderen 6 Probleme sind OFFEN."""
        result = millennium_problems_status()
        open_problems = [k for k, v in result['problems'].items()
                         if v['status'] == 'OFFEN']
        assert len(open_problems) == 6

    def test_summary_correct(self):
        """Zusammenfassung stimmt."""
        result = millennium_problems_status()
        summary = result['summary']
        assert summary['total'] == 7
        assert summary['solved'] == 1
        assert summary['open'] == 6

    def test_each_problem_has_prize(self):
        """Jedes Problem hat ein Preisgeld von 1.000.000 USD."""
        result = millennium_problems_status()
        for name, prob in result['problems'].items():
            assert prob['prize_usd'] == 1_000_000, f"{name}: falsches Preisgeld"

    def test_each_problem_structure(self):
        """Jedes Problem hat alle Pflichtfelder."""
        result = millennium_problems_status()
        required_fields = ['name', 'year_posed', 'field', 'status', 'prize_usd', 'description']
        for name, prob in result['problems'].items():
            for field in required_fields:
                assert field in prob, f"{name}: Feld '{field}' fehlt"

    def test_total_prize_remaining(self):
        """Verbleibendes Preisgeld = 6 × 1.000.000 = 6.000.000."""
        result = millennium_problems_status()
        assert result['summary']['total_prize_remaining_usd'] == 6_000_000


# =============================================================================
# EDGE-CASE TESTS
# =============================================================================

class TestEdgeCases:
    """Edge-Case Tests für ungewöhnliche Eingaben."""

    def test_goldbach_n4(self):
        """Kleinster gültiger Fall: 4 = 2 + 2."""
        result = goldbach_verification_range(4)
        assert result['all_verified'] is True

    def test_tsp_single_city(self):
        """TSP mit einer Stadt."""
        result = traveling_salesman_greedy([[0]])
        assert result['tour'] == [0]
        assert result['length'] == 0.0

    def test_tsp_two_cities(self):
        """TSP mit zwei Städten."""
        result = traveling_salesman_greedy([[0, 5], [5, 0]])
        assert len(result['tour']) == 2
        assert result['length'] == 10.0  # Hin- und Rückweg

    def test_complexity_n_large(self):
        """complexity_comparison für n=100 (2^100 und 100! → inf)."""
        result = complexity_comparison(100)
        # 2^100 ist sehr groß aber endlich repräsentierbar
        assert result['complexities']['O(2^n)'] == 2 ** 100
        # 100! ist > 20!, also inf
        assert result['complexities']['O(n!)'] == float('inf')

    def test_yang_mills_1d_field(self):
        """Yang-Mills für 1D-Feld."""
        F = np.array([1.0, 2.0, 3.0])
        action = yang_mills_action(F)
        # S = (1/4) * (1+4+9) = 14/4 = 3.5
        assert abs(action - 3.5) < 1e-12

    def test_navier_stokes_minimal_grid(self):
        """Navier-Stokes mit minimalem 5×5-Gitter."""
        result = navier_stokes_2d_simple(5, 5, dt=0.001, T=0.005, nu=0.1)
        assert isinstance(result, dict)
        assert result['u'].shape == (5, 5)

    def test_riemann_zeros_single(self):
        """Genau 1 Nullstelle berechnen."""
        zeros = riemann_zeros_mpmath(1)
        assert len(zeros) == 1
        assert zeros[0]['n'] == 1
