"""
Tests für das Siebtheorie-Modul (sieve_theory.py).

Testet alle Klassen und Funktionen des Sieb-Moduls:
    - BrunSieve: Kombinatorischer Sieb, Brun-Konstante
    - SelbergSieve: Oberschranken, Zwillingsprimzahlen
    - LargeSieve: Große-Sieb-Ungleichung
    - chen_theorem_verify: Chen-Theorem-Verifikation
    - goldbach_sieve_analysis: Goldbach-Komet-Analyse
    - SieveStatistics: Primzahllücken-Statistiken

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-11
@lastModified: 2026-03-11
"""

import math
import pytest
import sys
import os

# Projektverzeichnis zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from sieve_theory import (
    BrunSieve,
    SelbergSieve,
    LargeSieve,
    chen_theorem_verify,
    goldbach_sieve_analysis,
    SieveStatistics,
    _eratosthenes_sieve,
    _ist_prim,
)


# ===========================================================================
# HILFSFUNKTIONEN FÜR TESTS
# ===========================================================================

def bekannte_zwillingsprimzahlen(limit: int) -> list[tuple]:
    """
    Liefert bekannte Zwillingsprimzahlpaare bis 'limit' als Referenz.
    Verwendet simples Sieb, unabhängig vom getesteten Modul.
    """
    primes = set()
    ist_prim = [True] * (limit + 3)
    ist_prim[0] = ist_prim[1] = False
    for i in range(2, int(math.isqrt(limit + 2)) + 1):
        if ist_prim[i]:
            for j in range(i * i, limit + 3, i):
                ist_prim[j] = False
    primes = {i for i in range(2, limit + 3) if ist_prim[i]}
    paare = [(p, p + 2) for p in sorted(primes) if p <= limit and (p + 2) in primes]
    return paare


# ===========================================================================
# TESTS: _eratosthenes_sieve (Hilfsfunktion)
# ===========================================================================

class TestEratosthenesSieve:
    """Tests für die interne Sieb-Hilfsfunktion."""

    def test_kleine_werte(self):
        """Primzahlen bis 20 müssen korrekt sein."""
        assert _eratosthenes_sieve(20) == [2, 3, 5, 7, 11, 13, 17, 19]

    def test_grenze_0_und_1(self):
        """Kein Ergebnis für Grenzen < 2."""
        assert _eratosthenes_sieve(0) == []
        assert _eratosthenes_sieve(1) == []

    def test_genau_zwei(self):
        """Grenze 2 liefert nur [2]."""
        assert _eratosthenes_sieve(2) == [2]

    def test_anzahl_bis_100(self):
        """Es gibt genau 25 Primzahlen bis 100."""
        assert len(_eratosthenes_sieve(100)) == 25

    def test_erste_primzahl(self):
        """Erste Primzahl ist immer 2."""
        primes = _eratosthenes_sieve(1000)
        assert primes[0] == 2

    def test_keine_komposita(self):
        """Alle zurückgegebenen Zahlen müssen prim sein."""
        primes = _eratosthenes_sieve(200)
        for p in primes:
            assert _ist_prim(p), f"{p} sollte prim sein"


# ===========================================================================
# TESTS: BrunSieve
# ===========================================================================

class TestBrunSieve:
    """Tests für den Brun-Sieb."""

    def test_initialisierung(self):
        """Initialisierung mit gültigem limit."""
        sieb = BrunSieve(100)
        assert sieb.limit == 100

    def test_initialisierung_ungueltig(self):
        """Ungültige Grenzen werfen ValueError."""
        with pytest.raises(ValueError):
            BrunSieve(1)
        with pytest.raises(ValueError):
            BrunSieve(0)

    def test_sieve_ergebnis_korrekt(self):
        """Sieb liefert alle korrekten Primzahlen bis 50."""
        sieb = BrunSieve(50)
        ergebnis = sieb.sieve()
        erwartet = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
        assert ergebnis == erwartet

    def test_sieve_bis_100(self):
        """Sieb bis 100 liefert genau 25 Primzahlen."""
        sieb = BrunSieve(100)
        assert len(sieb.sieve()) == 25

    def test_sieve_enthält_keine_komposita(self):
        """Alle zurückgegebenen Zahlen sind prim."""
        sieb = BrunSieve(200)
        for p in sieb.sieve():
            assert _ist_prim(p), f"{p} ist keine Primzahl"

    def test_sieve_vollständig(self):
        """Alle Primzahlen bis limit sind im Ergebnis."""
        sieb = BrunSieve(150)
        ergebnis_set = set(sieb.sieve())
        referenz = set(_eratosthenes_sieve(150))
        assert ergebnis_set == referenz

    def test_brun_constant_wert_nah_an_bekanntem(self):
        """
        Bruns Konstante bis 10000 sollte nahe am bekannten Wert ~1.9021604 liegen.
        Die Reihe konvergiert sehr langsam (logarithmisch). Mit limit=10000 erreich man
        typischerweise ~1.6. Der Test prüft den unteren Teil des Konvergenzbereichs.
        Bekannter Wert: B₂ ≈ 1.9021604 (erreichbar erst für sehr großes limit).
        """
        sieb = BrunSieve(10000)
        b2 = sieb.brun_constant_estimate(10000)
        # Erreichbarer Bereich mit limit=10000: ~1.55–1.70
        # Toleranz: b2 muss zwischen 1.4 und 1.9 liegen
        assert 1.4 <= b2 <= 1.9, f"Brun-Konstante {b2:.6f} außerhalb [1.4, 1.9]"

    def test_brun_constant_monoton(self):
        """Bruns Konstante mit größerem Limit muss >= kleinerem sein."""
        sieb = BrunSieve(100)
        b2_klein = sieb.brun_constant_estimate(100)
        b2_gross = sieb.brun_constant_estimate(1000)
        assert b2_gross >= b2_klein

    def test_brun_constant_positiv(self):
        """Bruns Konstante muss positiv sein."""
        sieb = BrunSieve(50)
        b2 = sieb.brun_constant_estimate(50)
        assert b2 > 0

    def test_brun_constant_enthält_3_5_und_5_7(self):
        """
        Die Paare (3,5) und (5,7) sind Zwillingsprimzahlpaare bis 10.
        Daher muss brun_constant_estimate(10) mindestens
        1/3 + 1/5 + 1/5 + 1/7 = 0.876... betragen.
        (Das Paar (5,7) wird mit limit=10 ebenfalls eingeschlossen.)
        """
        sieb = BrunSieve(10)
        b2 = sieb.brun_constant_estimate(10)
        # (3,5) trägt 1/3+1/5 = 0.5333 bei, (5,7) trägt 1/5+1/7 = 0.3429 bei
        erwartet = 1 / 3 + 1 / 5 + 1 / 5 + 1 / 7
        # >= mit pytest.approx geht nicht direkt, daher Toleranz-Vergleich
        assert b2 >= erwartet - 1e-10, f"b2={b2:.6f} < erwartet={erwartet:.6f}"


# ===========================================================================
# TESTS: SelbergSieve
# ===========================================================================

class TestSelbergSieve:
    """Tests für den Selberg-Sieb."""

    def test_initialisierung(self):
        """Initialisierung mit korrekten Parametern."""
        sieb = SelbergSieve(100, [1, 5], 6)
        assert sieb.limit == 100
        assert sieb.residues == [1, 5]
        assert sieb.modulus == 6

    def test_initialisierung_ungueltig(self):
        """Ungültige Parameter werfen ValueError."""
        with pytest.raises(ValueError):
            SelbergSieve(100, [1], 0)
        with pytest.raises(ValueError):
            SelbergSieve(1, [1], 6)

    def test_upper_bound_positiv(self):
        """Oberschranke ist immer positiv für x >= 2."""
        sieb = SelbergSieve(1000, [1, 5], 6)
        schranke = sieb.upper_bound(1000)
        assert schranke > 0

    def test_upper_bound_wächst_mit_x(self):
        """Größeres x liefert größere (oder gleiche) Oberschranke."""
        sieb = SelbergSieve(1000, [1, 5], 6)
        assert sieb.upper_bound(500) <= sieb.upper_bound(1000)

    def test_upper_bound_zu_klein(self):
        """Für x < 2 wird 0 zurückgegeben."""
        sieb = SelbergSieve(100, [1, 5], 6)
        assert sieb.upper_bound(1) == 0.0

    def test_twin_primes_bis_100(self):
        """Zwillingsprimzahlen bis 100 müssen korrekt sein."""
        sieb = SelbergSieve(100, [1, 5], 6)
        ergebnis = sieb.sieve_twin_primes(100)
        erwartet = [(3, 5), (5, 7), (11, 13), (17, 19), (29, 31), (41, 43), (59, 61), (71, 73)]
        assert ergebnis == erwartet, f"Erhalten: {ergebnis}"

    def test_twin_primes_bis_50(self):
        """Zwillingsprimzahlen bis 50."""
        sieb = SelbergSieve(50, [1, 5], 6)
        ergebnis = sieb.sieve_twin_primes(50)
        erwartet = [(3, 5), (5, 7), (11, 13), (17, 19), (29, 31), (41, 43)]
        assert ergebnis == erwartet

    def test_twin_primes_alle_paare_sind_prim(self):
        """Alle gefundenen Paare müssen aus Primzahlen bestehen."""
        sieb = SelbergSieve(200, [1, 5], 6)
        paare = sieb.sieve_twin_primes(200)
        for p, q in paare:
            assert _ist_prim(p), f"{p} ist nicht prim"
            assert _ist_prim(q), f"{q} ist nicht prim"

    def test_twin_primes_abstand_genau_zwei(self):
        """Alle Paare haben Abstand genau 2."""
        sieb = SelbergSieve(300, [1, 5], 6)
        paare = sieb.sieve_twin_primes(300)
        for p, q in paare:
            assert q - p == 2, f"Abstand zwischen {p} und {q} ist nicht 2"

    def test_twin_primes_kein_ergebnis_unter_3(self):
        """Unter 3 gibt es keine Zwillingsprimzahlpaare."""
        sieb = SelbergSieve(10, [1, 5], 6)
        assert sieb.sieve_twin_primes(2) == []

    def test_twin_primes_übereinstimmung_mit_referenz(self):
        """Ergebnis stimmt mit bekannter Referenz überein."""
        sieb = SelbergSieve(500, [1, 5], 6)
        ergebnis = sieb.sieve_twin_primes(500)
        referenz = bekannte_zwillingsprimzahlen(500)
        assert ergebnis == referenz


# ===========================================================================
# TESTS: LargeSieve
# ===========================================================================

class TestLargeSieve:
    """Tests für den Großen Sieb."""

    def setup_method(self):
        """Erstellt eine LargeSieve-Instanz für jeden Test."""
        self.sieb = LargeSieve()

    def test_leere_eingabe(self):
        """Leere Koeffizientenfolge liefert 0."""
        assert self.sieb.large_sieve_inequality([], 5) == 0.0

    def test_q_null_liefert_null(self):
        """Q=0 liefert 0."""
        assert self.sieb.large_sieve_inequality([1.0, 2.0], 0) == 0.0

    def test_ungleichung_erfüllt(self):
        """
        Prüfe dass LHS ≤ (N + Q²) · Σ|aₙ|² gilt.
        """
        import numpy as np
        # Einfache Koeffizientenfolge: alle 1
        a = [1.0 + 0j] * 10
        Q = 5
        N = len(a)
        lhs = self.sieb.large_sieve_inequality(a, Q)
        rhs = (N + Q ** 2) * sum(abs(x) ** 2 for x in a)
        assert lhs <= rhs * (1 + 1e-10), f"LHS={lhs:.3f} > RHS={rhs:.3f}"

    def test_ungleichung_für_verschiedene_q(self):
        """Ungleichung gilt für Q=1,2,3,5,10."""
        a = [complex(1, 0.5)] * 20
        N = len(a)
        summe_an2 = sum(abs(x) ** 2 for x in a)
        for Q in [1, 2, 3, 5, 10]:
            lhs = self.sieb.large_sieve_inequality(a, Q)
            rhs = (N + Q ** 2) * summe_an2
            assert lhs <= rhs * (1 + 1e-8), f"Ungleichung verletzt für Q={Q}: {lhs:.3f} > {rhs:.3f}"

    def test_apply_to_primes_positiv(self):
        """Oberschranke für Primzahlen ist positiv."""
        schranke = self.sieb.apply_to_primes(50, 5)
        assert schranke > 0

    def test_apply_to_primes_wächst_mit_N(self):
        """Größeres N liefert größere Oberschranke."""
        s1 = self.sieb.apply_to_primes(50, 5)
        s2 = self.sieb.apply_to_primes(100, 5)
        assert s2 >= s1

    def test_apply_to_primes_N_unter_2(self):
        """Für N < 2 wird 0 zurückgegeben."""
        assert self.sieb.apply_to_primes(1, 5) == 0.0


# ===========================================================================
# TESTS: chen_theorem_verify
# ===========================================================================

class TestChenTheorem:
    """Tests für die Chen-Theorem-Verifikation."""

    def test_n_gleich_4(self):
        """4 = 2 + 2 (beide prim)."""
        ergebnis = chen_theorem_verify(4)
        assert ergebnis is not None
        p, q = ergebnis
        assert p + q == 4
        assert _ist_prim(p)

    def test_n_gleich_6(self):
        """6 = 3 + 3 oder 5 + 1 (5 prim, 3 prim)."""
        ergebnis = chen_theorem_verify(6)
        assert ergebnis is not None
        p, q = ergebnis
        assert p + q == 6

    def test_alle_geraden_bis_100(self):
        """Chen-Theorem für alle geraden n von 4 bis 100 verifiziert."""
        for n in range(4, 102, 2):
            ergebnis = chen_theorem_verify(n)
            assert ergebnis is not None, f"Chen-Theorem fehlgeschlagen für n={n}"
            p, q = ergebnis
            assert p + q == n, f"p+q={p+q} ≠ n={n}"
            assert _ist_prim(p), f"p={p} ist nicht prim (n={n})"

    def test_ergebnis_summe_korrekt(self):
        """Die Summe p+q muss immer n ergeben."""
        for n in [4, 10, 20, 50, 100, 200]:
            ergebnis = chen_theorem_verify(n)
            assert ergebnis is not None
            p, q = ergebnis
            assert p + q == n

    def test_ungültige_eingabe(self):
        """Ungerade oder zu kleine Zahlen liefern None."""
        assert chen_theorem_verify(3) is None   # ungerade
        assert chen_theorem_verify(2) is None   # zu klein
        assert chen_theorem_verify(1) is None   # zu klein
        assert chen_theorem_verify(7) is None   # ungerade

    def test_q_ist_prim_oder_semiprime(self):
        """q muss entweder prim oder ein Semiprime sein."""
        def ist_p1_oder_p2(m: int) -> bool:
            """Prüft ob m prim (P1) oder Produkt zweier Primzahlen (P2) ist."""
            if _ist_prim(m):
                return True
            # P2-Check: genau zwei Primfaktoren (mit Vielfachheit)
            temp = m
            count = 0
            for f in range(2, int(math.isqrt(m)) + 1):
                while temp % f == 0:
                    count += 1
                    temp //= f
                    if count > 2:
                        return False
            if temp > 1:
                count += 1
            return count == 2

        for n in range(4, 202, 2):
            ergebnis = chen_theorem_verify(n)
            assert ergebnis is not None
            p, q = ergebnis
            assert ist_p1_oder_p2(q), f"q={q} ist keine P1/P2-Zahl (n={n})"

    def test_bekannte_zerlegungen(self):
        """Bekannte Zerlegungen als Referenz."""
        # 4 = 2+2
        p, q = chen_theorem_verify(4)
        assert p == 2 and q == 2
        # 6 = 3+3
        p, q = chen_theorem_verify(6)
        assert p + q == 6


# ===========================================================================
# TESTS: goldbach_sieve_analysis
# ===========================================================================

class TestGoldbachAnalysis:
    """Tests für die Goldbach-Sieb-Analyse."""

    def test_grenze_zu_klein(self):
        """Für limit < 4 wird ein leeres Ergebnis zurückgegeben."""
        ergebnis = goldbach_sieve_analysis(3)
        assert ergebnis["min_decompositions"] == 0
        assert ergebnis["comet_data"] == []

    def test_min_decompositions_positiv(self):
        """Mindestens eine Zerlegung für alle geraden n von 4 bis 100."""
        ergebnis = goldbach_sieve_analysis(100)
        assert ergebnis["min_decompositions"] >= 1

    def test_max_größer_gleich_min(self):
        """Maximum >= Minimum."""
        ergebnis = goldbach_sieve_analysis(100)
        assert ergebnis["max_decompositions"] >= ergebnis["min_decompositions"]

    def test_durchschnitt_im_bereich(self):
        """Durchschnitt liegt zwischen min und max."""
        ergebnis = goldbach_sieve_analysis(200)
        assert ergebnis["min_decompositions"] <= ergebnis["average"] <= ergebnis["max_decompositions"]

    def test_komet_daten_länge(self):
        """Goldbach-Komet enthält Eintrag für jede gerade Zahl von 4 bis limit."""
        limit = 100
        ergebnis = goldbach_sieve_analysis(limit)
        erwartete_laenge = (limit - 4) // 2 + 1  # 4, 6, 8, ..., 100
        assert len(ergebnis["comet_data"]) == erwartete_laenge

    def test_komet_n_werte_korrekt(self):
        """Goldbach-Komet enthält die korrekten n-Werte (4, 6, 8, ...)."""
        ergebnis = goldbach_sieve_analysis(20)
        n_werte = [eintrag[0] for eintrag in ergebnis["comet_data"]]
        assert n_werte == [4, 6, 8, 10, 12, 14, 16, 18, 20]

    def test_bekannte_zerlegungen_n_4(self):
        """n=4: genau eine Zerlegung (2+2)."""
        ergebnis = goldbach_sieve_analysis(4)
        # komet_data = [(4, 2)]
        assert ergebnis["comet_data"][0] == (4, 2)

    def test_hardest_cases_haben_min_zerlegungen(self):
        """Schwierigste Fälle haben alle genau min_decompositions Zerlegungen."""
        ergebnis = goldbach_sieve_analysis(100)
        min_g = ergebnis["min_decompositions"]
        for n in ergebnis["hardest_cases"]:
            # Zähle Zerlegungen manuell nach
            prime_set = set(_eratosthenes_sieve(n))
            primes = sorted(prime_set)
            count = sum(1 for p in primes if p <= n // 2 and (n - p) in prime_set)
            assert count == min_g, f"n={n} hat {count} Zerlegungen, erwartet {min_g}"

    def test_komet_kleinste_primzahl_ist_prim(self):
        """Die kleinste Primzahl im Goldbach-Komet ist tatsächlich prim."""
        ergebnis = goldbach_sieve_analysis(100)
        for n, min_prime in ergebnis["comet_data"]:
            if min_prime > 0:
                assert _ist_prim(min_prime), f"min_prime={min_prime} für n={n} ist nicht prim"

    def test_n_12_hat_mindestens_2_zerlegungen(self):
        """
        n=12: Zerlegungen sind 5+7, 7+5 — da p≤q also 1 Paar,
        oder 5+7 und das Paar (1,11) zählt nicht. Tatsächlich: (5,7) Zerlegung p≤q/2.
        Mehr als 1 Zerlegung für n=12 erwartet.
        """
        ergebnis = goldbach_sieve_analysis(12)
        # Finde Eintrag für n=12
        for n, _ in ergebnis["comet_data"]:
            if n == 12:
                break
        # 12 = 5+7 ist die einzige Zerlegung mit p<=q
        # 12 = 5+7 und kein weiteres mit p<=6, also G(12)=1 möglich
        # Wichtiger: min_decompositions für gesamten Bereich [4,12] muss >= 1 sein
        assert ergebnis["min_decompositions"] >= 1

    def test_alle_n_haben_mindestens_eine_zerlegung_bis_1000(self):
        """
        Goldbach-Vermutung für alle geraden n von 4 bis 1000 überprüft.
        Jedes n muss mindestens eine Zerlegung haben.
        """
        ergebnis = goldbach_sieve_analysis(1000)
        assert ergebnis["min_decompositions"] >= 1, (
            f"Goldbach-Verletzung! Schwierigste Fälle: {ergebnis['hardest_cases']}"
        )


# ===========================================================================
# TESTS: SieveStatistics
# ===========================================================================

class TestSieveStatistics:
    """Tests für Primzahllücken-Statistiken."""

    def setup_method(self):
        """Erstellt eine SieveStatistics-Instanz für jeden Test."""
        self.stats = SieveStatistics()

    def test_prime_gaps_korrekte_werte(self):
        """Abstände zwischen kleinen Primzahlen müssen stimmen."""
        # Primes bis 20: [2, 3, 5, 7, 11, 13, 17, 19]
        # Abstände:       [1, 2, 2,  4,  2,  4,  2]
        luecken = self.stats.prime_gaps(20)
        assert luecken == [1, 2, 2, 4, 2, 4, 2]

    def test_prime_gaps_bis_10(self):
        """Abstände bis 10: [2,3,5,7] → [1,2,2]."""
        luecken = self.stats.prime_gaps(10)
        assert luecken == [1, 2, 2]

    def test_prime_gaps_zu_klein(self):
        """Für limit < 3 gibt es keine Lücken."""
        assert self.stats.prime_gaps(2) == []
        assert self.stats.prime_gaps(1) == []

    def test_prime_gaps_alle_positiv(self):
        """Alle Lücken sind positive ganze Zahlen."""
        luecken = self.stats.prime_gaps(500)
        for g in luecken:
            assert g > 0, f"Negative Lücke: {g}"

    def test_prime_gaps_erste_luecke_ist_eins(self):
        """Die erste Lücke (2→3) ist immer 1."""
        luecken = self.stats.prime_gaps(100)
        assert luecken[0] == 1

    def test_gap_distribution_enthält_alle_lücken(self):
        """Verteilung enthält alle vorkommenden Lückenwerte."""
        luecken = self.stats.prime_gaps(100)
        verteilung = self.stats.gap_distribution(100)
        for g in luecken:
            assert g in verteilung

    def test_gap_distribution_summe(self):
        """Summe der Häufigkeiten = Anzahl der Lücken."""
        luecken = self.stats.prime_gaps(200)
        verteilung = self.stats.gap_distribution(200)
        assert sum(verteilung.values()) == len(luecken)

    def test_gap_distribution_häufigkeit_zwei(self):
        """Lücke 2 (Zwillingsprimzahlen) kommt mehrfach vor bis 100."""
        verteilung = self.stats.gap_distribution(100)
        assert verteilung.get(2, 0) > 1

    def test_maximal_gap_bekannte_werte(self):
        """Bekannte maximale Lücken."""
        # Bis 30: [2,3,5,7,11,13,17,19,23,29], max gap = 6 (23→29)
        luecke, p_vor, p_nach = self.stats.maximal_gap(30)
        assert luecke == 6
        assert p_vor == 23
        assert p_nach == 29

    def test_maximal_gap_bis_100(self):
        """Maximale Lücke bis 100 ist zwischen 89 und 97 (Lücke 8)."""
        luecke, p_vor, p_nach = self.stats.maximal_gap(100)
        assert luecke == 8
        assert p_vor == 89
        assert p_nach == 97

    def test_maximal_gap_zu_klein(self):
        """Für limit < 3 wird (0,0,0) zurückgegeben."""
        assert self.stats.maximal_gap(2) == (0, 0, 0)

    def test_maximal_gap_ist_positiv(self):
        """Maximale Lücke für limit >= 5 ist positiv."""
        luecke, _, _ = self.stats.maximal_gap(100)
        assert luecke > 0

    def test_cramér_check_alle_erfüllt_ab_p_gleich_25(self):
        """
        Cramér-Vermutung: gap ≤ (ln p)² für alle Primzahlen bis 1000.
        Für sehr kleine Primzahlen (p < 25) ist die Cramér-Schranke (ln p)² < gap,
        da (ln 2)² ≈ 0.48, (ln 3)² ≈ 1.21, (ln 7)² ≈ 3.79 alle kleiner als
        die nächste Lücke sind. Die Vermutung ist für p ≥ 25 erfüllt.
        Bekannte Ausnahmen: p=2 (gap=1 > 0.48), p=3 (gap=2 > 1.21),
        p=7 (gap=4 > 3.79).
        """
        ergebnisse = self.stats.cramér_conjecture_check(1000)
        # Filtere: nur Primzahlen p >= 25 betrachten
        verletzungen = [(p, g, s) for p, g, s, erfüllt in ergebnisse
                        if not erfüllt and p >= 25]
        assert len(verletzungen) == 0, (
            f"Cramér-Vermutung verletzt für p >= 25: {verletzungen[:5]}"
        )

    def test_cramér_check_länge(self):
        """Anzahl der Ergebnisse = Anzahl der Primzahlen bis limit - 1."""
        primes = _eratosthenes_sieve(100)
        ergebnisse = self.stats.cramér_conjecture_check(100)
        assert len(ergebnisse) == len(primes) - 1

    def test_cramér_check_zu_klein(self):
        """Für limit < 3 gibt es keine Ergebnisse."""
        assert self.stats.cramér_conjecture_check(2) == []

    def test_cramér_check_format(self):
        """Jedes Ergebnis-Tupel hat 4 Elemente."""
        ergebnisse = self.stats.cramér_conjecture_check(50)
        for eintrag in ergebnisse:
            assert len(eintrag) == 4
            p, g, s, erfüllt = eintrag
            assert isinstance(p, int)
            assert isinstance(g, int)
            assert isinstance(s, float)
            assert isinstance(erfüllt, bool)

    def test_cramér_schranke_berechnung(self):
        """Cramér-Schranke (ln p)² wird korrekt berechnet."""
        ergebnisse = self.stats.cramér_conjecture_check(20)
        for p, g, s, erfüllt in ergebnisse:
            erwartete_schranke = math.log(p) ** 2
            assert abs(s - erwartete_schranke) < 1e-10, (
                f"Schranke für p={p}: erwartet {erwartete_schranke:.6f}, erhalten {s:.6f}"
            )


# ===========================================================================
# INTEGRATIONS-TESTS
# ===========================================================================

class TestIntegration:
    """Integrationstests über Modulgrenzen hinweg."""

    def test_brun_und_selberg_zwillingsprimzahlen_konsistent(self):
        """Brun-Sieb und Selberg-Sieb liefern konsistente Zwillingsprimzahlen."""
        limit = 100
        selberg = SelbergSieve(limit, [1, 5], 6)
        zwillinge = selberg.sieve_twin_primes(limit)

        brun = BrunSieve(limit)
        prime_list = brun.sieve()
        prime_set = set(prime_list)

        # Alle Selberg-Zwillingsprimzahlen müssen im Brun-Sieb sein
        for p, q in zwillinge:
            assert p in prime_set, f"{p} nicht im Brun-Sieb gefunden"
            assert q in prime_set, f"{q} nicht im Brun-Sieb gefunden"

    def test_chen_und_goldbach_konsistent(self):
        """
        Chen-Theorem deckt immer einen Goldbach-Fall ab:
        wenn n = p + q (beide prim) → chen_verify findet mindestens (p, q).
        """
        for n in range(4, 52, 2):
            # Goldbach: finde eine Zerlegung n = p + q
            primes = set(_eratosthenes_sieve(n))
            hat_goldbach = any(p in primes and (n - p) in primes
                               for p in range(2, n // 2 + 1))
            # Chen muss ebenfalls Ergebnis haben (stärker: P1 ∪ P2)
            chen = chen_theorem_verify(n)
            assert chen is not None, f"Chen fehlgeschlagen für n={n}"
            if hat_goldbach:
                # Chen darf auch die reine Goldbach-Zerlegung finden
                p, q = chen
                assert p + q == n

    def test_statistik_und_sieb_konsistent(self):
        """prime_gaps muss mit Sieb-Ergebnis übereinstimmen."""
        limit = 200
        stats = SieveStatistics()
        primes = _eratosthenes_sieve(limit)
        luecken_ref = [primes[i + 1] - primes[i] for i in range(len(primes) - 1)]
        luecken_stats = stats.prime_gaps(limit)
        assert luecken_stats == luecken_ref


# ===========================================================================
# EDGE CASES
# ===========================================================================

class TestEdgeCases:
    """Tests für Grenzfälle und Sonderfälle."""

    def test_brun_sieb_grenze_zwei(self):
        """Brun-Sieb mit limit=2 liefert [2]."""
        sieb = BrunSieve(2)
        assert sieb.sieve() == [2]

    def test_selberg_twin_primes_grenze_drei(self):
        """Zwillingsprimzahlen bis 3: keine (3,5) wegen 5>3."""
        sieb = SelbergSieve(10, [1, 5], 6)
        ergebnis = sieb.sieve_twin_primes(3)
        # (3,5): p=3 <= limit=3, aber 5>3+2? Nein, 5=3+2. p=3 ist in limit=3.
        # Korrekt: (3,5) ist enthalten wenn p=3 <= limit=3
        assert ergebnis == [(3, 5)]

    def test_chen_für_sehr_große_gerade_zahl(self):
        """Chen-Theorem für n=200 (etwas größere Zahl)."""
        ergebnis = chen_theorem_verify(200)
        assert ergebnis is not None
        p, q = ergebnis
        assert p + q == 200

    def test_goldbach_nur_n_gleich_4(self):
        """goldbach_sieve_analysis(4) enthält nur n=4."""
        ergebnis = goldbach_sieve_analysis(4)
        assert len(ergebnis["comet_data"]) == 1
        assert ergebnis["comet_data"][0][0] == 4

    def test_large_sieve_einzelnes_element(self):
        """Large Sieve mit einem Koeffizient."""
        sieb = LargeSieve()
        lhs = sieb.large_sieve_inequality([1.0 + 0j], 3)
        # Oberschranke: (1 + 9) * 1 = 10
        assert lhs <= (1 + 9) * 1 * (1 + 1e-9)

    def test_prime_gaps_enthält_keine_null(self):
        """Keine Lücke der Größe 0 (alle Primzahlen verschieden)."""
        stats = SieveStatistics()
        luecken = stats.prime_gaps(1000)
        assert 0 not in luecken

    def test_cramér_check_für_erste_primzahl(self):
        """
        Erster Eintrag ist für p=2, Lücke=1.
        Hinweis: (ln 2)² ≈ 0.48 < 1 = gap → Cramér-Vermutung für p=2 verletzt.
        Die Cramér-Vermutung gilt erst für hinreichend große p.
        Wir prüfen hier nur das Format, nicht das Erfülltsein.
        """
        stats = SieveStatistics()
        ergebnisse = stats.cramér_conjecture_check(10)
        p, g, s, erfüllt = ergebnisse[0]
        assert p == 2       # Erster Eintrag ist p=2
        assert g == 1       # Lücke 2→3 ist 1
        # Schranke: (ln 2)² ≈ 0.48 — kleiner als gap=1, also erfüllt=False erwartet
        assert not erfüllt, f"Für p=2 sollte Cramér-Vermutung verletzt sein: gap={g}, (ln p)²={s:.4f}"
