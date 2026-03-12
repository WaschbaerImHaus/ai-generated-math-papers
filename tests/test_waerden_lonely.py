"""
Tests für van_der_waerden.py und lonely_runner.py.

Umfassende Test-Suite mit ≥80 Tests für:
  - VanDerWaerdenNumbers: Bekannte Werte, AP-Erkennung, Schranken
  - ArithmeticProgressionColoring: SAT-Kodierung, Erfüllbarkeit
  - LonelyRunnerConjecture: Verifikation für n=2..6
  - LonelyRunnerVerifier: Exakte Zeitpunktsuche

@author: Michael Fuhrmann
@date: 2026-03-12
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import pytest
import math
from fractions import Fraction
from typing import List

# Module importieren
from van_der_waerden import VanDerWaerdenNumbers, ArithmeticProgressionColoring
from lonely_runner import LonelyRunnerConjecture, LonelyRunnerVerifier


# ============================================================
# Hilfsfunktionen für Tests
# ============================================================

def make_coloring(assignment: dict, N: int) -> List[int]:
    """
    Erzeugt eine Färbungsliste aus einem Dictionary {i: Farbe}.

    @param assignment: {1-basierter Index: Farbe (0-basiert)}
    @param N: Größe des Bereichs
    @return: Färbungsliste (0-basiert, Länge N)
    """
    return [assignment.get(i + 1, 0) for i in range(N)]


# ============================================================
# Klasse: TestArithmeticProgressionBasics
# ============================================================

class TestArithmeticProgressionBasics:
    """Tests für grundlegende AP-Erkennung und -Prüfung."""

    def setup_method(self):
        self.vdw = VanDerWaerdenNumbers()

    # --- is_arithmetic_progression ---

    def test_ap_length2_is_ap(self):
        """Zwei Elemente bilden immer eine AP."""
        assert self.vdw.is_arithmetic_progression([1, 3]) is True

    def test_ap_length1_is_ap(self):
        """Ein einzelnes Element ist trivial eine AP."""
        assert self.vdw.is_arithmetic_progression([5]) is True

    def test_ap_classic_135(self):
        """[1,3,5] ist AP mit Differenz 2."""
        assert self.vdw.is_arithmetic_progression([1, 3, 5]) is True

    def test_ap_classic_246(self):
        """[2,4,6] ist AP mit Differenz 2."""
        assert self.vdw.is_arithmetic_progression([2, 4, 6]) is True

    def test_ap_classic_1_2_3(self):
        """[1,2,3] ist AP mit Differenz 1."""
        assert self.vdw.is_arithmetic_progression([1, 2, 3]) is True

    def test_not_ap_124(self):
        """[1,2,4] ist KEINE AP (Differenzen: 1,2)."""
        assert self.vdw.is_arithmetic_progression([1, 2, 4]) is False

    def test_not_ap_135_plus_9(self):
        """[1,3,5,9] ist KEINE AP (letzte Differenz 4≠2)."""
        assert self.vdw.is_arithmetic_progression([1, 3, 5, 9]) is False

    def test_ap_constant_all_equal(self):
        """[3,3,3] ist AP mit Differenz 0."""
        assert self.vdw.is_arithmetic_progression([3, 3, 3]) is True

    def test_ap_decreasing(self):
        """[9,6,3] ist AP mit negativer Differenz -3."""
        assert self.vdw.is_arithmetic_progression([9, 6, 3]) is True

    def test_ap_length4(self):
        """[1,4,7,10] ist AP mit Differenz 3."""
        assert self.vdw.is_arithmetic_progression([1, 4, 7, 10]) is True

    def test_not_ap_random(self):
        """[1,5,6] ist keine AP."""
        assert self.vdw.is_arithmetic_progression([1, 5, 6]) is False

    # --- detect_ap ---

    def test_detect_ap_finds_135(self):
        """Erkennt AP [1,3,5] in Sequenz [1,2,3,4,5]."""
        aps = self.vdw.detect_ap([1, 2, 3, 4, 5], k=3)
        assert [1, 3, 5] in aps

    def test_detect_ap_empty_in_random(self):
        """In [1,2,4] gibt es keine 3-term AP."""
        aps = self.vdw.detect_ap([1, 2, 4], k=3)
        # Keine AP mit Länge 3 und positivem d
        for ap in aps:
            assert not self.vdw.is_arithmetic_progression(ap) or len(ap) < 3 or False
        # Kurze Liste kann [1,2], [2,4], [1,4] enthalten — aber keine 3-term AP
        for ap in aps:
            assert len(ap) == 3

    def test_detect_ap_multiple(self):
        """[1,2,3,4,5] enthält mehrere 3-term APs."""
        aps = self.vdw.detect_ap([1, 2, 3, 4, 5], k=3)
        assert len(aps) >= 3  # mind. [1,2,3], [1,3,5], [2,3,4], [3,4,5], [2,4,6?], ...

    def test_detect_ap_length4(self):
        """Erkennt 4-term AP [1,4,7,10] in Sequenz."""
        aps = self.vdw.detect_ap([1, 2, 4, 7, 10], k=4)
        assert [1, 4, 7, 10] in aps

    # --- is_ap_free_coloring ---

    def test_ap_free_coloring_trivial_true(self):
        """[0,1] (2-Färbung, nur 2 Elemente) ist trivial AP-frei für k=3."""
        assert self.vdw.is_ap_free_coloring([0, 1], k=3) is True

    def test_ap_free_coloring_false_example(self):
        """Färbung [0,0,0] enthält 3-term AP {1,2,3}."""
        assert self.vdw.is_ap_free_coloring([0, 0, 0], k=3) is False

    def test_ap_free_coloring_correct_8(self):
        """Bekannte AP-freie 2-Färbung von {1,...,8} für k=3."""
        # Aus SAT-Lösung: W(3;2)=9, also {1,...,8} ist 2-färbbar ohne 3-AP
        # Beispiel: [1,1,0,0,1,1,0,0] (0-basiert, gefunden via SAT)
        coloring = [1, 1, 0, 0, 1, 1, 0, 0]
        assert self.vdw.is_ap_free_coloring(coloring, k=3) is True

    def test_ap_contains_monochromatic(self):
        """[0,1,0,1,0] — prüfe ob 3-term AP vorhanden."""
        coloring = [0, 1, 0, 1, 0]
        # {0,2,4} = Indizes 0,2,4 → Farbe 0,0,0 → monochromatisch
        result = self.vdw.is_ap_free_coloring(coloring, k=3)
        assert result is False


# ============================================================
# Klasse: TestVanDerWaerdenKnownValues
# ============================================================

class TestVanDerWaerdenKnownValues:
    """Tests für bekannte Van-der-Waerden-Zahlen."""

    def setup_method(self):
        self.vdw = VanDerWaerdenNumbers()

    def test_lookup_w32_equals_9(self):
        """W(3;2) = 9: Kleinste N sodass jede 2-Färbung eine mono 3-AP enthält."""
        assert self.vdw.lookup(3, 2) == 9

    def test_lookup_w42_equals_35(self):
        """W(4;2) = 35 (bekannter Wert)."""
        assert self.vdw.lookup(4, 2) == 35

    def test_lookup_w52_equals_178(self):
        """W(5;2) = 178 (bekannter Wert)."""
        assert self.vdw.lookup(5, 2) == 178

    def test_lookup_w33_equals_27(self):
        """W(3;3) = 27 (bekannter Wert)."""
        assert self.vdw.lookup(3, 3) == 27

    def test_lookup_w34_equals_293(self):
        """W(3;4) = 293 (bekannter Wert)."""
        assert self.vdw.lookup(3, 4) == 293

    def test_lookup_unknown_returns_none(self):
        """W(4;4) ist unbekannt → lookup gibt None zurück."""
        assert self.vdw.lookup(4, 4) is None

    def test_lookup_w22_none(self):
        """W(2;2) ist nicht in der Tabelle."""
        result = self.vdw.lookup(2, 2)
        assert result is None


# ============================================================
# Klasse: TestVanDerWaerdenSATCompute
# ============================================================

class TestVanDerWaerdenSATCompute:
    """Tests für SAT-basierte Berechnung von W(k;r)."""

    def setup_method(self):
        self.vdw = VanDerWaerdenNumbers()

    def test_compute_w32_via_sat(self):
        """W(3;2) = 9 durch SAT-Berechnung."""
        result = self.vdw.compute(3, 2, max_N=15)
        assert result == 9

    def test_compute_w33_via_sat(self):
        """W(3;3) = 27 durch SAT-Berechnung."""
        result = self.vdw.compute(3, 3, max_N=30)
        assert result == 27

    def test_compute_w22_via_sat(self):
        """W(2;2) = 3: Jede 2-Färbung von {1,2,3} hat 2-term AP (Schubfach)."""
        # W(2;2) = 3: {1,2} kann AP-frei gefärbt werden (1→0, 2→1).
        # {1,2,3} mit 2 Farben: mind. 2 Elemente gleiche Farbe = 2-AP.
        result = self.vdw.compute(2, 2, max_N=10)
        assert result == 3

    def test_lower_bound_w32(self):
        """W(3;2) ≥ 9: {1,...,8} ist 2-färbbar ohne 3-term AP."""
        assert self.vdw.verify_lower_bound(9, k=3, r=2) is True

    def test_upper_bound_w32(self):
        """W(3;2) ≤ 9: {1,...,9} hat keine 2-Färbung ohne 3-term AP."""
        assert self.vdw.verify_upper_bound(9, k=3, r=2) is True

    def test_lower_bound_w33(self):
        """W(3;3) ≥ 27: {1,...,26} ist 3-färbbar ohne 3-term AP."""
        assert self.vdw.verify_lower_bound(27, k=3, r=3) is True

    def test_sat_coloring_8_exists(self):
        """Für N=8, k=3, r=2 gibt es eine gültige AP-freie Färbung."""
        coloring = self.vdw.find_ap_free_coloring(8, k=3, r=2)
        assert coloring is not None
        assert len(coloring) == 8

    def test_sat_coloring_9_not_exists(self):
        """Für N=9, k=3, r=2 gibt es KEINE AP-freie 2-Färbung."""
        coloring = self.vdw.find_ap_free_coloring(9, k=3, r=2)
        assert coloring is None

    def test_sat_coloring_is_valid(self):
        """Zurückgegebene AP-freie Färbung ist tatsächlich gültig."""
        coloring = self.vdw.find_ap_free_coloring(8, k=3, r=2)
        assert coloring is not None
        assert self.vdw.is_ap_free_coloring(coloring, k=3) is True

    def test_sat_coloring_26_exists_k3_r3(self):
        """W(3;3)=27: {1,...,26} hat eine gültige 3-Färbung ohne 3-AP."""
        coloring = self.vdw.find_ap_free_coloring(26, k=3, r=3)
        assert coloring is not None

    def test_compute_uses_cache(self):
        """Zweifacher Aufruf gibt identisches Ergebnis (Cache korrekt)."""
        r1 = self.vdw.compute(3, 2, max_N=15)
        r2 = self.vdw.compute(3, 2, max_N=15)
        assert r1 == r2 == 9

    def test_compute_returns_none_if_not_found(self):
        """compute() gibt None zurück falls max_N zu klein für unbekannte Parameter."""
        # W(4;4) ist unbekannt, max_N=10 zu klein
        result = self.vdw.compute(4, 4, max_N=10)
        assert result is None


# ============================================================
# Klasse: TestBerlekampBound
# ============================================================

class TestBerlekampBound:
    """Tests für die Berlekamp-Untergrenze."""

    def setup_method(self):
        self.vdw = VanDerWaerdenNumbers()

    def test_berlekamp_p2(self):
        """Berlekamp p=2: W(3;2) > 2·2^2 = 8 → W(3;2) ≥ 9 ✓."""
        lb = self.vdw.berlekamp_lower_bound(2)
        assert lb == 8
        # W(3;2) = 9 > 8 ✓
        assert self.vdw.lookup(3, 2) == 9
        assert 9 > lb

    def test_berlekamp_p3(self):
        """Berlekamp p=3: W(4;2) > 3·2^3 = 24."""
        lb = self.vdw.berlekamp_lower_bound(3)
        assert lb == 24

    def test_berlekamp_p5(self):
        """Berlekamp p=5: W(6;2) > 5·2^5 = 160."""
        lb = self.vdw.berlekamp_lower_bound(5)
        assert lb == 160

    def test_berlekamp_p7(self):
        """Berlekamp p=7: W(8;2) > 7·2^7 = 896."""
        lb = self.vdw.berlekamp_lower_bound(7)
        assert lb == 896

    def test_berlekamp_formula(self):
        """Berlekamp-Formel: lb(p) = p * 2^p."""
        for p, expected in [(2, 8), (3, 24), (5, 160), (7, 896), (11, 11 * 2048)]:
            assert self.vdw.berlekamp_lower_bound(p) == expected

    def test_berlekamp_p2_consistent_with_w32(self):
        """W(3;2) = 9 > Berlekamp(2) = 8 (konsistent)."""
        # Berlekamp: W(p+1;2) > p*2^p für Primzahl p
        # p=2: W(3;2) > 8 → W(3;2) ≥ 9 ✓
        assert self.vdw.lookup(3, 2) == 9
        result = self.vdw.compute(3, 2, max_N=12)
        assert result == 9
        assert result > self.vdw.berlekamp_lower_bound(2)


# ============================================================
# Klasse: TestArithmeticProgressionColoringDirect
# ============================================================

class TestArithmeticProgressionColoringDirect:
    """Direkte Tests der ArithmeticProgressionColoring-Klasse."""

    def test_var_id_basic(self):
        """Variablen-IDs sind korrekt berechnet."""
        coder = ArithmeticProgressionColoring(N=5, k=3, r=2)
        # var(i,c) = (i-1)*r + c
        assert coder._var(1, 1) == 1
        assert coder._var(1, 2) == 2
        assert coder._var(2, 1) == 3
        assert coder._var(2, 2) == 4

    def test_compute_aps_k3_n5(self):
        """Für N=5, k=3 werden alle 3-term APs korrekt gefunden."""
        coder = ArithmeticProgressionColoring(N=5, k=3, r=2)
        aps = coder._aps
        # Erwartete APs: [1,2,3],[1,3,5],[2,3,4],[2,4,6?],[3,4,5]...
        # Nur APs mit letztem Element ≤ 5
        for ap in aps:
            assert len(ap) == 3
            assert ap[-1] <= 5
            assert ap[1] - ap[0] == ap[2] - ap[1]  # konstante Differenz

    def test_compute_aps_correct_count(self):
        """Anzahl der APs ist korrekt."""
        coder = ArithmeticProgressionColoring(N=6, k=3, r=2)
        # Alle 3-term APs in {1,...,6}:
        # d=1: (1,2,3),(2,3,4),(3,4,5),(4,5,6) → 4
        # d=2: (1,3,5),(2,4,6) → 2
        assert len(coder._aps) == 6

    def test_cnf_satisfiable_for_n8_k3_r2(self):
        """N=8, k=3, r=2 ist SAT-erfüllbar."""
        coder = ArithmeticProgressionColoring(N=8, k=3, r=2)
        result = coder.solve()
        assert result is not None

    def test_cnf_unsatisfiable_for_n9_k3_r2(self):
        """N=9, k=3, r=2 ist UNSAT."""
        coder = ArithmeticProgressionColoring(N=9, k=3, r=2)
        result = coder.solve()
        assert result is None

    def test_extract_coloring_dict_valid(self):
        """extract_coloring_dict gibt gültiges Dict zurück."""
        coder = ArithmeticProgressionColoring(N=8, k=3, r=2)
        d = coder.extract_coloring_dict()
        assert d is not None
        assert len(d) == 8
        assert all(k in range(1, 9) for k in d.keys())
        assert all(v in (1, 2) for v in d.values())

    def test_extract_coloring_dict_none_for_unsat(self):
        """extract_coloring_dict gibt None für UNSAT zurück."""
        coder = ArithmeticProgressionColoring(N=9, k=3, r=2)
        d = coder.extract_coloring_dict()
        assert d is None

    def test_solved_coloring_is_ap_free(self):
        """Gelöste Färbung besteht tatsächlich den AP-Freiheits-Test."""
        coder = ArithmeticProgressionColoring(N=8, k=3, r=2)
        coloring = coder.solve()
        assert coloring is not None
        vdw = VanDerWaerdenNumbers()
        # Konvertiere: Farben 1,2 → 0,1 für is_ap_free_coloring
        converted = [c - 1 for c in coloring]
        assert vdw.is_ap_free_coloring(converted, k=3) is True


# ============================================================
# Klasse: TestLonelyRunnerConjecture
# ============================================================

class TestLonelyRunnerConjecture:
    """Tests für die Lonely-Runner-Vermutung."""

    # --- Grundlegende Konstruktion ---

    def test_init_basic(self):
        """Grundlegende Initialisierung funktioniert."""
        lrc = LonelyRunnerConjecture([0, 1, 2])
        assert lrc.n == 3

    def test_init_normalization(self):
        """Normalisierung auf v[0]=0 funktioniert korrekt."""
        lrc = LonelyRunnerConjecture([3, 5, 7])
        # Normiert: 0, 2, 4
        assert lrc.velocities[0] == 0

    def test_init_duplicate_velocities_raises(self):
        """Doppelte Geschwindigkeiten werfen ValueError."""
        with pytest.raises(ValueError):
            LonelyRunnerConjecture([0, 1, 1, 2])

    def test_threshold_n2(self):
        """Schwellwert für n=2 ist 1/2."""
        lrc = LonelyRunnerConjecture([0, 1])
        assert lrc.threshold == Fraction(1, 2)

    def test_threshold_n4(self):
        """Schwellwert für n=4 ist 1/4."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3])
        assert lrc.threshold == Fraction(1, 4)

    def test_relative_velocities(self):
        """Relative Geschwindigkeiten (ohne 0) korrekt."""
        lrc = LonelyRunnerConjecture([0, 1, 3])
        assert 0 not in lrc.relative_velocities

    # --- circle_distance ---

    def test_circle_distance_zero(self):
        """Abstand 0 auf dem Kreis."""
        lrc = LonelyRunnerConjecture([0, 1])
        assert abs(lrc.circle_distance(0.0)) < 1e-9

    def test_circle_distance_half(self):
        """Abstand 1/2 auf dem Kreis (gegenüberliegende Seite)."""
        lrc = LonelyRunnerConjecture([0, 1])
        assert abs(lrc.circle_distance(0.5) - 0.5) < 1e-9

    def test_circle_distance_quarter(self):
        """Abstand 1/4."""
        lrc = LonelyRunnerConjecture([0, 1])
        assert abs(lrc.circle_distance(0.25) - 0.25) < 1e-9

    def test_circle_distance_three_quarters(self):
        """3/4 ergibt Abstand 1/4 (kürzere Seite)."""
        lrc = LonelyRunnerConjecture([0, 1])
        assert abs(lrc.circle_distance(0.75) - 0.25) < 1e-9

    def test_circle_distance_fraction_exact(self):
        """Exakte Fraction-Abstandsberechnung."""
        lrc = LonelyRunnerConjecture([0, 1])
        d = lrc.circle_distance_fraction(Fraction(1, 3))
        assert d == Fraction(1, 3)

    def test_circle_distance_fraction_two_thirds(self):
        """2/3 ergibt Abstand 1/3."""
        lrc = LonelyRunnerConjecture([0, 1])
        d = lrc.circle_distance_fraction(Fraction(2, 3))
        assert d == Fraction(1, 3)

    # --- n=2: Trivialfall ---

    def test_n2_runner0_lonely_at_half(self):
        """n=2, v=[0,1]: Läufer 0 ist lonely bei t=1/2."""
        lrc = LonelyRunnerConjecture([0, 1])
        # Bei t=1/2: Position von Läufer 1 ist 1/2, Abstand = 1/2 ≥ 1/2 ✓
        assert lrc.is_lonely_at_time(0, 0.5) is True

    def test_n2_runner0_lonely_at_half_v2(self):
        """n=2, v=[0,2]: Läufer 0 ist lonely bei t=1/4 (relative Pos = 1/2)."""
        lrc = LonelyRunnerConjecture([0, 2])
        # Bei t=1/4: relative Pos = 2*1/4 = 1/2, Abstand = 1/2 ≥ 1/2 ✓
        assert lrc.is_lonely_at_time(0, 0.25) is True

    def test_n2_not_lonely_at_zero(self):
        """Bei t=0 sind alle Läufer an Position 0 — nicht lonely."""
        lrc = LonelyRunnerConjecture([0, 1])
        # t=0: Abstand = 0 < 1/2
        assert lrc.is_lonely_at_time(0, 0.0) is False

    def test_n2_verifier_trivial(self):
        """LonelyRunnerVerifier.verify_n2 gibt immer True zurück."""
        verifier = LonelyRunnerVerifier()
        ok, t = verifier.verify_n2(v1=1)
        assert ok is True
        assert t == Fraction(1, 2)

    def test_n2_verifier_v3(self):
        """verify_n2 für v1=3: t=1/6."""
        verifier = LonelyRunnerVerifier()
        ok, t = verifier.verify_n2(v1=3)
        assert ok is True
        assert t == Fraction(1, 6)

    # --- n=3 ---

    def test_n3_conjecture_verified(self):
        """n=3 Vermutung: Alle Läufer haben einen Loneliness-Zeitpunkt."""
        lrc = LonelyRunnerConjecture([0, 1, 2])
        # Für Läufer 0: Suche t mit ∥t·1∥ ≥ 1/3 und ∥t·2∥ ≥ 1/3
        # t=1/3: ∥1/3∥ = 1/3 ≥ 1/3 ✓, ∥2/3∥ = 1/3 ≥ 1/3 ✓
        assert lrc.is_lonely_at_time_exact(0, Fraction(1, 3)) is True

    def test_n3_simulation_finds_time(self):
        """Simulation findet Loneliness-Zeitpunkt für n=3."""
        lrc = LonelyRunnerConjecture([0, 1, 3])
        t = lrc.find_lonely_time_simulation(0, t_max=10.0, steps=10000)
        assert t is not None
        assert t > 0

    def test_n3_min_distance(self):
        """Minimaler Abstand von Läufer 0 bei bekanntem guten t."""
        lrc = LonelyRunnerConjecture([0, 1, 2])
        # t=1/3: beide anderen Läufer bei 1/3 und 2/3 → Abstand je 1/3
        d = lrc.min_distance_at_time(0, 1.0 / 3.0)
        assert d >= 1.0 / 3.0 - 1e-9

    # --- n=4 ---

    def test_n4_conjecture_runner0(self):
        """n=4: Läufer 0 hat Loneliness-Zeitpunkt."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3])
        t = lrc.find_lonely_time_simulation(0, t_max=5.0, steps=50000)
        assert t is not None

    def test_n4_is_lonely_exact(self):
        """n=4, v=[0,1,2,3]: Prüfe konkrete Zeitpunkte."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3])
        # Versuche t=1/4
        # ∥1/4·1∥=1/4 ≥ 1/4 ✓, ∥1/4·2∥=∥1/2∥=1/2≥1/4 ✓, ∥1/4·3∥=∥3/4∥=1/4≥1/4 ✓
        assert lrc.is_lonely_at_time_exact(0, Fraction(1, 4)) is True

    def test_n4_blocki_sketch_not_empty(self):
        """blocki_sketch_n4 gibt nicht-leeren String zurück."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3])
        sketch = lrc.blocki_sketch_n4()
        assert len(sketch) > 100
        assert "n=4" in sketch or "Läufer" in sketch

    # --- n=5 ---

    def test_n5_conjecture_verified_example(self):
        """n=5: Bekanntes Beispiel [0,1,2,3,4] — Läufer 0 lonely bei t=1/5."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3, 4])
        # t=1/5: ∥1/5∥=1/5≥1/5✓, ∥2/5∥=2/5≥1/5✓, ∥3/5∥=2/5≥1/5✓, ∥4/5∥=1/5≥1/5✓
        assert lrc.is_lonely_at_time_exact(0, Fraction(1, 5)) is True

    def test_n5_threshold_correct(self):
        """n=5: Schwellwert = 1/5."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3, 4])
        assert lrc.threshold == Fraction(1, 5)

    # --- n=6 ---

    def test_n6_threshold(self):
        """n=6: Schwellwert = 1/6."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3, 4, 5])
        assert lrc.threshold == Fraction(1, 6)

    def test_n6_conjecture_simulation(self):
        """n=6: Simulation findet Loneliness-Zeitpunkt für Läufer 0."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 4, 5, 6])
        t = lrc.find_lonely_time_simulation(0, t_max=20.0, steps=100000)
        assert t is not None

    # --- Diophantische Reformulierung ---

    def test_diophantine_description_not_empty(self):
        """Diophantische Reformulierung gibt Text zurück."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3])
        desc = lrc.describe_diophantine_reformulation()
        assert "1/n" in desc or "Conjecture" in desc.lower() or "n =" in desc

    def test_diophantine_contains_threshold(self):
        """Beschreibung enthält Schwellwert."""
        lrc = LonelyRunnerConjecture([0, 1, 2, 3])
        desc = lrc.describe_diophantine_reformulation()
        assert "1/4" in desc or "0.25" in desc


# ============================================================
# Klasse: TestLonelyRunnerVerifier
# ============================================================

class TestLonelyRunnerVerifier:
    """Tests für LonelyRunnerVerifier."""

    def setup_method(self):
        self.verifier = LonelyRunnerVerifier()

    def test_verify_n2_v1_result(self):
        """verify_n2 für v1=1 gibt (True, 1/2)."""
        ok, t = self.verifier.verify_n2(v1=1)
        assert ok is True
        assert t == Fraction(1, 2)

    def test_verify_n2_v2_result(self):
        """verify_n2 für v1=2 gibt (True, 1/4)."""
        ok, t = self.verifier.verify_n2(v1=2)
        assert ok is True
        assert t == Fraction(1, 4)

    def test_verify_runner_exact_n3_simple(self):
        """Exakte Verifikation für n=3, v=[0,1,2], Läufer 0."""
        ok, t = self.verifier.verify_runner_exact([0, 1, 2], runner_idx=0)
        assert ok is True
        assert t is not None
        # Prüfe dass t tatsächlich ein Loneliness-Zeitpunkt ist
        lrc = LonelyRunnerConjecture([0, 1, 2])
        assert lrc.is_lonely_at_time_exact(0, t) is True

    def test_verify_runner_exact_n4(self):
        """Exakte Verifikation für n=4, v=[0,1,2,3]."""
        ok, t = self.verifier.verify_runner_exact([0, 1, 2, 3], runner_idx=0)
        assert ok is True

    def test_verify_all_runners_n3(self):
        """Alle Läufer für n=3 verifiziert."""
        results = self.verifier.verify_all_runners([0, 1, 2])
        # Alle 3 Läufer sollten verifiziert sein
        for i, (ok, t) in results.items():
            assert ok is True, f"Läufer {i} nicht verifiziert"

    def test_generate_test_cases_n2(self):
        """Test-Case-Generierung für n=2."""
        cases = self.verifier.generate_test_cases(n=2, max_vel=5)
        assert len(cases) > 0
        for case in cases:
            assert len(case) == 2
            assert len(set(case)) == 2  # paarweise verschieden

    def test_generate_test_cases_n3(self):
        """Test-Case-Generierung für n=3."""
        cases = self.verifier.generate_test_cases(n=3, max_vel=5)
        assert len(cases) > 0
        for case in cases:
            assert len(case) == 3


# ============================================================
# Klasse: TestEdgeCases
# ============================================================

class TestEdgeCases:
    """Edge-Cases und Sonderfälle."""

    def setup_method(self):
        self.vdw = VanDerWaerdenNumbers()

    def test_ap_free_single_element(self):
        """Einzelnes Element ist trivial AP-frei."""
        assert self.vdw.is_ap_free_coloring([0], k=3) is True

    def test_ap_free_two_elements(self):
        """Zwei Elemente sind AP-frei für k=3."""
        assert self.vdw.is_ap_free_coloring([0, 1], k=3) is True

    def test_sat_k2_trivial(self):
        """k=2: Jede Färbung mit 2 gleichen Elementen hat 2-AP."""
        # N=2, k=2, r=1: Alle Elemente haben Farbe 1 → 2-AP {1,2}
        coder = ArithmeticProgressionColoring(N=2, k=2, r=1)
        result = coder.solve()
        assert result is None  # UNSAT: 2 Elemente gleiche Farbe = 2-AP

    def test_lonely_runner_large_velocities(self):
        """Große Geschwindigkeiten: n=2, v=[0,100]."""
        lrc = LonelyRunnerConjecture([0, 100])
        # t = 1/200: Relative Pos = 100*1/200 = 1/2 → Abstand = 1/2 ≥ 1/2 ✓
        assert lrc.is_lonely_at_time(0, 0.005) is True

    def test_lonely_runner_negative_velocities(self):
        """Negative Geschwindigkeiten werden korrekt normiert."""
        lrc = LonelyRunnerConjecture([-1, 0, 1])
        # Normiert auf [0, 1, 2]
        assert lrc.velocities[0] == 0
        assert 1 in lrc.velocities
        assert 2 in lrc.velocities

    def test_count_ap_free_colorings_small(self):
        """count_ap_free_colorings für kleine N."""
        # Für N=3, k=3, r=2: Alle 2^3=8 Färbungen prüfen
        # {0,0,0}, {1,1,1} haben AP → nicht gezählt
        # {0,0,1}: Indizes 0,1 gleich → aber k=3 braucht 3 gleiche
        count = self.vdw.count_ap_free_colorings(3, k=3, r=2)
        # {0,0,0} und {1,1,1} sind nicht AP-frei → 8-2=6 sind AP-frei
        assert count == 6

    def test_berlekamp_consistency(self):
        """Berlekamp-Untergrenze ist konsistent mit W(3;2)=9."""
        # W(p+1;2) > p*2^p für p=2: W(3;2) > 8 → W(3;2) ≥ 9 ✓
        lb = self.vdw.berlekamp_lower_bound(2)
        w32 = self.vdw.compute(3, 2, max_N=12)
        assert w32 is not None
        assert w32 > lb

    def test_gowers_bound_k3(self):
        """Gowers-Schranke für k=3 ist sehr groß."""
        bound = self.vdw.gowers_upper_bound(3)
        # 2^(2^12) = astronomisch groß
        assert bound > 10 ** 100

    def test_ap_detection_empty_sequence(self):
        """Leere Sequenz enthält keine AP."""
        aps = self.vdw.detect_ap([], k=3)
        assert aps == []

    def test_ap_detection_sequence_too_short(self):
        """Zu kurze Sequenz für k=4."""
        aps = self.vdw.detect_ap([1, 2, 3], k=4)
        assert aps == []

    def test_lonely_runner_verify_all_returns_dict(self):
        """verify_conjecture_for_all_runners gibt Dict zurück."""
        lrc = LonelyRunnerConjecture([0, 1, 2])
        result = lrc.verify_conjecture_for_all_runners(t_max=5.0, steps=5000)
        assert isinstance(result, dict)
        assert len(result) == 3
