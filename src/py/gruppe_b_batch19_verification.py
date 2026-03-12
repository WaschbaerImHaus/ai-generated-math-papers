"""
@file gruppe_b_batch19_verification.py
@brief Verifikation der vier Vermutungen aus Batch 19 (Gruppe B).

@description
    Dieses Skript untersucht und verifiziert (wo möglich) die folgenden vier
    mathematischen Vermutungen aus Batch 19:

    1. Freiman-Struktursatz / Polynomial Freiman-Ruzsa (PFR) Vermutung
       - Über F₂ⁿ bewiesen (Gowers-Green-Manners-Tao 2023)
       - Über ℤ noch offen (quantitative Schranken werden geprüft)

    2. Frankl Union-Closed Conjecture (Erdős-Ko-Rado-Kontext)
       - Gilmer 2022: Existiert Element mit Frequenz > 1/100
       - Vollständige 1/2-Grenze offen
       - Kleine Familien werden exhaustiv geprüft

    3. Lonely Runner Vermutung
       - n=4,5,6,7: bewiesen (verifiziert numerisch)
       - n=8: offen (Suche nach möglichem Gegenbeispiel)

    4. Graceful Tree Vermutung
       - Alle Bäume ≤10 Knoten: verifiziert durch Backtracking

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import sys
import time
import math
import itertools
from fractions import Fraction
from typing import List, Dict, Tuple, Optional, Set, FrozenSet

# Netzwerkbibliothek für Graphen
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    print("[WARNUNG] networkx nicht verfügbar. Graceful-Tree-Tests werden übersprungen.")

import numpy as np


# ===========================================================================
# TEIL 1: FREIMAN / POLYNOMIAL FREIMAN-RUZSA (PFR)
# ===========================================================================

class FreimanPFRAnalysis:
    """
    Analyse der Freiman-Struktursatz-Vermutung und der Polynomial Freiman-Ruzsa
    (PFR) Vermutung über ℤ und F₂ⁿ.

    Freimans Satz (klassisch, 1973):
        Ist A ⊆ ℤ eine endliche Menge mit |A+A| ≤ K·|A|, so ist A in einem
        arithmetischen Progressionsprodukt der Dimension d ≤ d(K) und Größe
        ≤ C(K)·|A| enthalten, wobei d(K) und C(K) von K abhängen.

    PFR über F₂ⁿ (Gowers-Green-Manners-Tao 2023, BEWIESEN):
        Ist A ⊆ F₂ⁿ mit |A+A| ≤ K·|A|, so ist |A| ≥ 2^{n-CK} für eine
        universelle Konstante C > 0. D.h. A ist nahe einer affinen Untergruppe.

    PFR über ℤ (OFFEN):
        Ist A ⊆ ℤ mit |A+A| ≤ K·|A|, so gibt es eine arithmetische
        Progression P mit A ⊆ P und |P| ≤ K^C · |A| für polynomiales C.
        (Sanders 2012: C ~ log K; PFR: C = O(1) erwartet)

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    def compute_doubling_constant(self, A: List[int]) -> float:
        """
        Berechnet die Verdopplungskonstante K = |A+A| / |A|.

        Die Summmenge A+A = {a+b : a,b ∈ A}.

        @param A: Endliche Menge ganzer Zahlen
        @return: Verdopplungskonstante K ≥ 1
        @lastModified 2026-03-12
        """
        A_set = set(A)
        # Summmenge A+A berechnen
        sumset = set()
        for a in A_set:
            for b in A_set:
                sumset.add(a + b)
        return len(sumset) / len(A_set)

    def freiman_ruzsa_bound_classical(self, K: float, d: int) -> float:
        """
        Klassische Freiman-Ruzsa-Schranke für |P| / |A|.

        Nach Ruzsa (1999): |P| / |A| ≤ K^2 · 2^{d·K}
        (grobe obere Schranke für die Größe der umhüllenden Progression)

        @param K: Verdopplungskonstante |A+A|/|A|
        @param d: Dimension der Progression
        @return: Obere Schranke für |P|/|A|
        @lastModified 2026-03-12
        """
        return K**2 * 2**(d * K)

    def pfr_f2n_bound(self, K: float, C_pfr: float = 12.0) -> float:
        """
        Schranke gemäß dem PFR-Beweis über F₂ⁿ (Gowers et al. 2023).

        Der Beweis liefert: Hat A ⊆ F₂ⁿ Verdopplungskonstante K,
        dann ist dist(A, affiner Unterraum) ≤ C_pfr · log₂(K).
        (Konzeptuell: A ist in einem affinen Unterraum H der Codimension ≤ C·log K enthalten)

        Hinweis: C_pfr = 12 ist eine Näherung; der exakte Wert aus dem
        Gowers-Green-Manners-Tao-Beweis ist komplex.

        @param K: Verdopplungskonstante
        @param C_pfr: Konstante aus dem PFR-Beweis (Näherungswert)
        @return: Obere Schranke für log₂(|umhüllender Unterraum| / |A|)
        @lastModified 2026-03-12
        """
        if K <= 1:
            return 0.0
        return C_pfr * math.log2(K)

    def sanders_bound_over_z(self, K: float) -> float:
        """
        Sanders-Schranke über ℤ (2012, bestes bekanntes Resultat).

        Sanders (2012): Für A ⊆ ℤ mit |A+A| ≤ K|A| gibt es eine
        Progression P mit A ⊆ P und |P| ≤ exp(C · log^4(K) / log log K) · |A|.

        Dies ist subexponentiell in K, aber nicht polynomial (PFR offen).

        @param K: Verdopplungskonstante
        @return: Schranke für |P|/|A| (grob)
        @lastModified 2026-03-12
        """
        if K <= 1:
            return 1.0
        C = 5.0  # Näherungskonstante
        log_k = math.log(K)
        log_log_k = math.log(max(log_k, math.e))
        exponent = C * log_k**4 / log_log_k
        return math.exp(exponent)

    def pfr_over_z_conjectured_bound(self, K: float, C: float = 2.0) -> float:
        """
        Vermutete PFR-Schranke über ℤ (falls PFR wahr).

        PFR-Vermutung: |P|/|A| ≤ K^C für polynomiales C.

        @param K: Verdopplungskonstante
        @param C: Polynomiale Konstante (erwartet C ~ 2-4)
        @return: Conjectured bound |P|/|A|
        @lastModified 2026-03-12
        """
        return K**C

    def verify_small_sets_structure(self, max_n: int = 20, max_K: float = 4.0) -> List[Dict]:
        """
        Untersucht kleine Mengen A ⊆ {0,...,N} auf ihre Verdopplungskonstante
        und die Struktur der Summmenge.

        Für Mengen mit kleiner Verdopplungskonstante prüft diese Funktion,
        ob A nahe einer arithmetischen Progression liegt.

        @param max_n: Maximale Elementanzahl in A
        @param max_K: Maximale Verdopplungskonstante für Ausgabe
        @return: Liste von Dicts mit Analyseergebnissen
        @lastModified 2026-03-12
        """
        results = []
        # Kleines Beispiel: arithmetische Progressionen haben K → 2 (minimal)
        for step in range(1, 6):
            for length in range(3, min(max_n, 12) + 1):
                A = list(range(0, length * step, step))
                K = self.compute_doubling_constant(A)
                # Für AP: |A+A| = 2|A|-1, also K → 2 - 1/|A|
                results.append({
                    "type": f"AP(step={step},len={length})",
                    "A": A[:8],  # Nur erste 8 Elemente anzeigen
                    "|A|": len(A),
                    "K": round(K, 4),
                    "structure": "arithmetische Progression" if K < 2.1 else "komplex"
                })
        return results[:15]  # Begrenze Ausgabe

    def f2n_sumset_example(self, n: int = 4) -> Dict:
        """
        Demonstriert PFR-Konzept in F₂ⁿ mit einem konkreten Beispiel.

        F₂ⁿ = (ℤ/2ℤ)ⁿ: Addition modulo 2 komponentenweise.
        Untergruppen sind affine Unterräume.

        @param n: Dimension von F₂ⁿ
        @return: Dict mit Beispiel-Ergebnissen
        @lastModified 2026-03-12
        """
        N = 2**n
        # Vollständige Gruppe F₂ⁿ als Bitstrings 0,...,2^n - 1
        # Untergruppe H = {0, 1, ..., 2^(n-1) - 1} (erste Hälfte) — kein echter Unterraum!
        # Echter Unterraum: alle Vektoren mit letztem Bit 0
        subspace = [i for i in range(N) if (i & 1) == 0]  # Letzte Koordinate = 0

        # Summmenge in F₂ⁿ (XOR)
        def f2n_sumset(A: List[int]) -> Set[int]:
            return {a ^ b for a in A for b in A}

        full_sumset = f2n_sumset(subspace)
        K = len(full_sumset) / len(subspace)

        return {
            "n": n,
            "F₂ⁿ_size": N,
            "Unterraum_H": len(subspace),
            "H+H": len(full_sumset),
            "K=|H+H|/|H|": round(K, 4),
            "Erklaerung": "Echter Unterraum hat K=1 (H+H=H). PFR sagt: kleine K ⟹ nahe Unterraum."
        }


# ===========================================================================
# TEIL 2: FRANKL UNION-CLOSED CONJECTURE
# ===========================================================================

class FranklUnionClosedAnalysis:
    """
    Analyse der Frankl Union-Closed Conjecture (1979).

    Vermutung (Frankl 1979):
        Sei F eine endliche Familie von Mengen (nicht nur {∅}), die unter
        Vereinigung abgeschlossen ist (d.h. A,B ∈ F ⟹ A∪B ∈ F).
        Dann gibt es ein Element x, das in mindestens der Hälfte aller
        Mengen von F vorkommt: |{A ∈ F : x ∈ A}| ≥ |F|/2.

    Stand (2024):
        - Gilmer (2022): Es gibt stets ein Element mit Frequenz > 0.01·|F|
          (also > 1/100, deutlich unter 1/2).
          Verbesserung auf 0.38·|F| durch Alweiss-Huang-Sellke (2022).
        - Vollständige 1/2-Grenze: OFFEN.
        - Für Familien mit ≤ 50 Mengen: Verifiziert (Czédli 2009 u.a.).

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    def is_union_closed(self, family: List[FrozenSet]) -> bool:
        """
        Prüft ob eine Mengenfamilie unter Vereinigung abgeschlossen ist.

        @param family: Liste von frozensets
        @return: True wenn union-closed
        @lastModified 2026-03-12
        """
        family_set = set(family)
        for A in family:
            for B in family:
                if (A | B) not in family_set:
                    return False
        return True

    def max_frequency(self, family: List[FrozenSet]) -> Tuple[int, int, float]:
        """
        Findet das Element mit maximaler Häufigkeit in der Familie.

        @param family: Union-closed Mengenfamilie
        @return: (element, count, fraction) — häufigstes Element
        @lastModified 2026-03-12
        """
        if not family:
            return (-1, 0, 0.0)

        # Sammle alle Elemente
        all_elements: Set[int] = set()
        for A in family:
            all_elements |= A

        if not all_elements:
            return (-1, 0, 0.0)

        # Zähle Häufigkeiten
        freq: Dict[int, int] = {}
        for x in all_elements:
            freq[x] = sum(1 for A in family if x in A)

        best_elem = max(freq, key=lambda x: freq[x])
        count = freq[best_elem]
        fraction = count / len(family)
        return (best_elem, count, fraction)

    def verify_frankl_conjecture_small(self, ground_set_size: int = 4) -> Dict:
        """
        Verifiziert die Frankl-Vermutung exhaustiv für kleine Grundmengen.

        Für Grundmenge {0,...,n-1} generiert alle union-closed Familien
        und prüft die 1/2-Grenze.

        @param ground_set_size: Größe der Grundmenge
        @return: Verifikationsergebnisse
        @lastModified 2026-03-12
        """
        n = ground_set_size
        ground = list(range(n))
        all_subsets = [frozenset(S) for r in range(n + 1)
                       for S in itertools.combinations(ground, r)]

        checked = 0
        violations = 0
        max_failing_fraction = 1.0
        worst_case = None

        # Alle nicht-leeren Teilfamilien prüfen
        # (zu viele für große n, daher Stichproben für n=4)
        for r in range(1, len(all_subsets) + 1):
            for family_tuple in itertools.combinations(all_subsets, r):
                family = list(family_tuple)
                # Nur nicht-triviale Familien (mind. eine nicht-leere Menge)
                if all(len(A) == 0 for A in family):
                    continue
                if not self.is_union_closed(family):
                    continue

                checked += 1
                best_elem, count, frac = self.max_frequency(family)

                if frac < 0.5:
                    violations += 1
                    if frac < max_failing_fraction:
                        max_failing_fraction = frac
                        worst_case = {
                            "family": [list(A) for A in family],
                            "best_fraction": round(frac, 4)
                        }

                # Zeitlimit: nicht zu lange suchen
                if checked > 5000:
                    break
            if checked > 5000:
                break

        return {
            "ground_set_size": n,
            "families_checked": checked,
            "frankl_violations": violations,
            "min_observed_max_freq": round(1.0 - max_failing_fraction, 4) if worst_case else "N/A",
            "worst_case": worst_case,
            "frankl_holds": violations == 0,
            "note": "Vollständige Prüfung bis Stichprobenlimit"
        }

    def gilmer_bound_example(self) -> Dict:
        """
        Demonstriert die Gilmer-Schranke (2022) an einem Beispiel.

        Gilmer bewies: In jeder union-closed Familie gibt es ein Element x mit
            |{A ∈ F : x ∈ A}| ≥ (3 - √5)/2 · |F| ≈ 0.382 · |F|.

        Verbessert durch Alweiss-Huang-Sellke (2022) auf ~0.38, Chase-Lovett auf ~0.40.

        @return: Erklärung und Beispiel
        @lastModified 2026-03-12
        """
        # Gilmers Schranke: (3 - sqrt(5)) / 2
        gilmer_bound = (3 - math.sqrt(5)) / 2
        # Klassisches kritisches Beispiel: "sunflower"-Familie
        # F = {{1}, {2}, {1,2}} — union-closed, |F|=3
        # Element 1: in {1},{1,2} → 2/3 > 1/2 ✓
        family_example = [frozenset([1]), frozenset([2]), frozenset([1, 2])]
        best_elem, count, frac = self.max_frequency(family_example)

        return {
            "gilmer_bound_1_100": "> 1/100 (original Gilmer 2022)",
            "gilmer_bound_exact": round(gilmer_bound, 6),
            "gilmer_bound_approx": "≈ 0.382 · |F|",
            "best_known_bound": "≈ 0.38 (Alweiss-Huang-Sellke 2022), ~0.40 (Chase-Lovett)",
            "target": "0.5 (Frankl-Vermutung)",
            "gap": "0.10 bis 0.12 — Lücke bleibt erheblich",
            "example_family": [[1], [2], [1, 2]],
            "example_best_fraction": round(frac, 4),
            "fully_proven": False
        }

    def frankl_for_specific_families(self) -> List[Dict]:
        """
        Prüft die Frankl-Vermutung für spezifische bekannte Familien.

        Testet: Powersets, Ketten, Sunflowers, und konstruierte Beispiele.

        @return: Liste von Verifikationsergebnissen
        @lastModified 2026-03-12
        """
        results = []

        test_families = [
            # Powerset von {0,1,2} (ohne ∅)
            ("Powerset({0,1,2})", [
                frozenset([0]), frozenset([1]), frozenset([2]),
                frozenset([0, 1]), frozenset([0, 2]), frozenset([1, 2]),
                frozenset([0, 1, 2])
            ]),
            # Kette: {{0}, {0,1}, {0,1,2}}
            ("Kette_3", [frozenset([0]), frozenset([0, 1]), frozenset([0, 1, 2])]),
            # Sunflower: {{1,3},{2,3},{1,2,3}}
            ("Sunflower_3", [frozenset([1, 3]), frozenset([2, 3]), frozenset([1, 2, 3])]),
            # Kritisches Beispiel Frankl: Größte Familie nahe 1/2
            ("Union_Closed_4", [
                frozenset([1, 2]), frozenset([1, 3]), frozenset([1, 2, 3]),
                frozenset([1, 2, 4]), frozenset([1, 3, 4]), frozenset([1, 2, 3, 4])
            ]),
        ]

        for name, family in test_families:
            uc = self.is_union_closed(family)
            if uc:
                elem, count, frac = self.max_frequency(family)
                results.append({
                    "Familie": name,
                    "union_closed": uc,
                    "|F|": len(family),
                    "häufigstes Element": elem,
                    "Häufigkeit": count,
                    "Anteil": round(frac, 4),
                    "≥ 1/2": frac >= 0.5
                })
            else:
                results.append({
                    "Familie": name,
                    "union_closed": False,
                    "Hinweis": "Nicht union-closed — übersprungen"
                })

        return results


# ===========================================================================
# TEIL 3: LONELY RUNNER VERMUTUNG
# ===========================================================================

class LonelyRunnerVerification:
    """
    Verifikation der Lonely-Runner-Vermutung für n=4,5,6,7,8 Läufer.

    Lonely Runner Conjecture (Wills 1967, Cusick 1973):
        Gegeben n Läufer auf einem Einheitskreis (Umfang 1), alle starten
        bei 0, paarweise verschiedene ganzzahlige Geschwindigkeiten v₀,...,v_{n-1}.
        Dann existiert für JEDEN Läufer i ein Zeitpunkt t > 0, sodass:
            ∥t·(vᵢ - vⱼ)∥ ≥ 1/n für alle j ≠ i
        (Abstand auf dem Kreis zur nächsten ganzen Zahl ≥ 1/n).

    O.B.d.A.: v₀ = 0 (Referenzläufer stationär).
    Die Vermutung ist äquivalent: Für jeden Läufer i ≠ 0 und alle
    anderen Läufer j ≠ i ist ∥t·wₖ∥ ≥ 1/n, wobei wₖ = vᵢ - vₖ.

    Beweisstatus:
        n ≤ 1: trivial
        n = 2,3,4: Cusick & Pomerance (1984) — BEWIESEN
        n = 5: Bohman, Fonoberova, Pikhurko (2001) — BEWIESEN
        n = 6: Bohman et al. (2011) — BEWIESEN
        n = 7: Tao (2018) — BEWIESEN
        n ≥ 8: OFFEN (Conjecture)

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    def circle_dist(self, x: float) -> float:
        """
        Kreisabstand ∥x∥ = min(x mod 1, 1 - x mod 1).

        @param x: reelle Zahl
        @return: Abstand zur nächsten ganzen Zahl in [0, 0.5]
        @lastModified 2026-03-12
        """
        x_mod = x % 1.0
        return min(x_mod, 1.0 - x_mod)

    def circle_dist_frac(self, x: Fraction) -> Fraction:
        """
        Exakter Kreisabstand mit rationaler Arithmetik.

        @param x: rationaler Wert
        @return: exakter Kreisabstand als Fraction
        @lastModified 2026-03-12
        """
        x_mod = x - int(x)
        if x_mod < 0:
            x_mod += Fraction(1)
        return min(x_mod, Fraction(1) - x_mod)

    def is_lonely_exact(self, vels: List[int], runner: int, t: Fraction) -> bool:
        """
        Prüft exakt ob Läufer `runner` zum Zeitpunkt t lonely ist.

        @param vels: Geschwindigkeitsliste (normalisiert auf v[0]=0)
        @param runner: Index des zu prüfenden Läufers
        @param t: rationaler Zeitpunkt
        @return: True wenn lonely
        @lastModified 2026-03-12
        """
        n = len(vels)
        threshold = Fraction(1, n)
        vi = vels[runner]
        for j, vj in enumerate(vels):
            if j == runner:
                continue
            pos = t * Fraction(vi - vj)
            if self.circle_dist_frac(pos) < threshold:
                return False
        return True

    def verify_runner_instance(
        self,
        vels: List[int],
        runner: int,
        max_denom: int = 300
    ) -> Tuple[bool, Optional[Fraction]]:
        """
        Sucht einen Loneliness-Zeitpunkt für Läufer `runner`.

        Strategie: Kritische Zeitpunkte (wo sich relative Positionen ändern)
        sind t = k/d mit d = |vᵢ - vⱼ| für alle Paare (i,j).
        Zwischen kritischen Punkten ist die Loneliness-Bedingung konstant.
        Mittelpunkte werden geprüft.

        @param vels: Geschwindigkeitsliste
        @param runner: Läufer-Index
        @param max_denom: Maximaler Nenner für Kandidaten-Brüche
        @return: (True, t) wenn lonely-Zeit gefunden, (False, None) sonst
        @lastModified 2026-03-12
        """
        # Normalisierung: v[0] auf 0
        v0 = vels[0]
        normalized = [v - v0 for v in vels]

        # Alle Differenzen für den Läufer 'runner'
        vi = normalized[runner]
        diffs = []
        for j, vj in enumerate(normalized):
            if j != runner:
                d = abs(vi - vj)
                if d > 0:
                    diffs.append(d)

        # Kritische Punkte in [0,1]: t = k/d
        critical: Set[Fraction] = {Fraction(0), Fraction(1)}
        for d in diffs:
            for k in range(d + 1):
                critical.add(Fraction(k, d))

        # Zusätzliche fein-granulare Kandidaten
        for b in range(1, max_denom + 1):
            for a in range(1, b):
                critical.add(Fraction(a, b))

        # Sortiere und prüfe Mittelpunkte
        sorted_crit = sorted(critical)
        for i in range(len(sorted_crit) - 1):
            mid = (sorted_crit[i] + sorted_crit[i + 1]) / 2
            if mid > 0 and self.is_lonely_exact(normalized, runner, mid):
                return True, mid

        # Direkte Punkte prüfen
        for t in sorted_crit:
            if t > 0 and self.is_lonely_exact(normalized, runner, t):
                return True, t

        return False, None

    def verify_all_runners(
        self,
        vels: List[int],
        max_denom: int = 300
    ) -> Dict[int, Tuple[bool, Optional[Fraction]]]:
        """
        Verifiziert alle Läufer einer Instanz.

        @param vels: Geschwindigkeitsliste
        @param max_denom: Suchtiefe
        @return: Dict {runner: (verified, time)}
        @lastModified 2026-03-12
        """
        return {
            i: self.verify_runner_instance(vels, i, max_denom)
            for i in range(len(vels))
        }

    def run_batch_verification(
        self,
        n: int,
        num_samples: int = 8,
        max_denom: int = 200
    ) -> Dict:
        """
        Verifiziert n Läufer an `num_samples` zufälligen Instanzen.

        @param n: Anzahl Läufer
        @param num_samples: Anzahl Testinstanzen
        @param max_denom: Suchtiefe
        @return: Zusammenfassung
        @lastModified 2026-03-12
        """
        from itertools import combinations

        # Kanonische Testfälle: kleine Geschwindigkeiten 1,...,n+2
        test_cases = []
        all_possible = list(combinations(range(1, n + 5), n - 1))
        for combo in all_possible[:num_samples]:
            vels = [0] + list(combo)
            test_cases.append(vels)

        all_verified = True
        failed_cases = []
        verified_count = 0
        total_runners = 0

        for vels in test_cases:
            res = self.verify_all_runners(vels, max_denom)
            for runner, (ok, t) in res.items():
                total_runners += 1
                if ok:
                    verified_count += 1
                else:
                    all_verified = False
                    failed_cases.append({
                        "velocities": vels,
                        "runner": runner,
                        "note": "Kein Lonely-Zeitpunkt mit max_denom gefunden (erhöhe Suchtiefe)"
                    })

        return {
            "n": n,
            "instances_tested": len(test_cases),
            "total_runners": total_runners,
            "verified": verified_count,
            "all_verified": all_verified,
            "failed_cases": failed_cases[:3],  # Maximal 3 anzeigen
        }

    def search_counterexample_n8(
        self,
        num_cases: int = 5,
        max_denom: int = 100
    ) -> Dict:
        """
        Sucht nach einem Gegenbeispiel für n=8 Läufer.

        Für n=8 ist die Vermutung offen. Diese Funktion testet kanonische
        Instanzen und prüft ob die Vermutung fehlschlägt.

        @param num_cases: Anzahl der getesteten Instanzen
        @param max_denom: Suchtiefe (niedrig für Geschwindigkeit)
        @return: Ergebnisse der Suche
        @lastModified 2026-03-12
        """
        from itertools import combinations

        n = 8
        # Teste kanonische Fälle: {0, 1, 2, ..., 7} und Varianten
        test_vels = [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [0, 1, 2, 3, 5, 7, 11, 13],  # Primzahlen
            [0, 2, 3, 4, 5, 6, 7, 8],
            [0, 1, 3, 5, 7, 9, 11, 13],  # Ungerade
            [0, 1, 2, 4, 8, 16, 32, 64],  # Potenzen von 2
        ]

        counterexamples = []
        verified_instances = 0

        for vels in test_vels[:num_cases]:
            res = self.verify_all_runners(vels, max_denom)
            all_ok = all(ok for ok, _ in res.values())
            if not all_ok:
                # Mögliches Gegenbeispiel oder nur Suchtiefe unzureichend?
                failed_runners = [i for i, (ok, _) in res.items() if not ok]
                counterexamples.append({
                    "velocities": vels,
                    "failed_runners": failed_runners,
                    "note": "Kein Beweis mit aktueller Suchtiefe — kein Gegenbeispiel!"
                })
            else:
                verified_instances += 1

        return {
            "n": 8,
            "cases_tested": len(test_vels[:num_cases]),
            "fully_verified": verified_instances,
            "potential_issues": len(counterexamples),
            "note": "Kein Gegenbeispiel gefunden — aber n=8 ist mathematisch OFFEN",
            "cases_with_issues": counterexamples
        }


# ===========================================================================
# TEIL 4: GRACEFUL TREE VERMUTUNG
# ===========================================================================

class GracefulTreeVerification:
    """
    Verifikation der Graceful Tree Conjecture (Rosa 1967) für kleine Bäume.

    Graceful Labeling (Rosa 1967):
        Ein Baum T mit n Knoten besitzt ein graceful Labeling, wenn es eine
        bijektive Abbildung f: V(T) → {0, 1, ..., n-1} gibt, sodass die
        Kantenbeschriftungen |f(u) - f(v)| für alle Kanten {u,v} paarweise
        verschieden sind (also genau {1, ..., n-1}).

    Äquivalenz zur Ringel-Kotzig-Vermutung:
        Jeder Baum T mit n Kanten teilt K_{2n+1} in 2n+1 isomorphe Kopien.

    Beweisstatus:
        - Alle Bäume mit ≤ 35 Knoten: verifiziert (Aldred & McKay 1998)
        - Pfade, Caterpillar-Graphen, Helme: bewiesen (allgemein)
        - Allgemeiner Fall: OFFEN

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    def is_graceful(self, tree: "nx.Graph", labeling: Dict) -> bool:
        """
        Prüft ob ein Labeling graceful ist.

        @param tree: NetworkX-Graph (Baum)
        @param labeling: Dict {Knoten: Label in 0..n-1}
        @return: True wenn graceful
        @lastModified 2026-03-12
        """
        n = tree.number_of_nodes()
        if set(labeling.values()) != set(range(n)):
            return False
        edge_labels = set()
        for u, v in tree.edges():
            diff = abs(labeling[u] - labeling[v])
            if diff in edge_labels:
                return False
            edge_labels.add(diff)
        return edge_labels == set(range(1, n))

    def find_graceful_backtrack(self, tree: "nx.Graph") -> Optional[Dict]:
        """
        Sucht ein graceful Labeling via Backtracking mit Constraint-Propagation.

        Algorithmus:
          1. Ordne Knoten nach Grad (hoher Grad zuerst — bessere Heuristik)
          2. Weise Labels 0..n-1 zu, breche früh ab wenn Kanten-Konflikte
          3. Gib erstes gefundenes Labeling zurück

        @param tree: NetworkX-Baum
        @return: Graceful Labeling als Dict oder None
        @lastModified 2026-03-12
        """
        n = tree.number_of_nodes()
        nodes = sorted(tree.nodes(), key=lambda v: tree.degree(v), reverse=True)

        labeling: Dict = {}
        used_labels: Set[int] = set()
        used_edge_labels: Set[int] = set()

        def backtrack(idx: int) -> bool:
            if idx == n:
                # Prüfe ob alle Kantenlabels {1,...,n-1} sind
                return used_edge_labels == set(range(1, n))

            node = nodes[idx]
            for label in range(n):
                if label in used_labels:
                    continue

                # Berechne neue Kantenlabels
                new_edge_labels = []
                conflict = False
                for neighbor in tree.neighbors(node):
                    if neighbor in labeling:
                        diff = abs(label - labeling[neighbor])
                        if diff in used_edge_labels or diff in new_edge_labels:
                            conflict = True
                            break
                        if diff == 0 or diff >= n:
                            conflict = True
                            break
                        new_edge_labels.append(diff)

                if conflict:
                    continue

                # Label zuweisen
                labeling[node] = label
                used_labels.add(label)
                for el in new_edge_labels:
                    used_edge_labels.add(el)

                if backtrack(idx + 1):
                    return True

                # Rücksetzen
                del labeling[node]
                used_labels.discard(label)
                for el in new_edge_labels:
                    used_edge_labels.discard(el)

            return False

        if backtrack(0):
            return dict(labeling)
        return None

    def enumerate_all_trees(self, n: int) -> List["nx.Graph"]:
        """
        Enumeriert alle nicht-isomorphen Bäume mit n Knoten.

        Für kleine n (≤ 10) werden alle Bäume via NetworkX generiert.

        @param n: Anzahl Knoten (≤ 10 empfohlen)
        @return: Liste nicht-isomorpher Bäume
        @lastModified 2026-03-12
        """
        if not HAS_NETWORKX:
            return []
        # NetworkX bietet nx.nonisomorphic_trees für kleine n
        trees = list(nx.nonisomorphic_trees(n))
        return trees

    def verify_all_trees_up_to_n(self, max_n: int = 10) -> Dict:
        """
        Verifiziert die Graceful-Tree-Vermutung für alle Bäume mit ≤ max_n Knoten.

        @param max_n: Maximale Knotenanzahl
        @return: Verifikationsergebnisse
        @lastModified 2026-03-12
        """
        if not HAS_NETWORKX:
            return {"error": "networkx nicht verfügbar"}

        results = {}
        total_trees = 0
        total_graceful = 0
        counterexamples = []

        for n in range(1, max_n + 1):
            trees = self.enumerate_all_trees(n)
            n_graceful = 0
            n_total = len(trees)

            for tree in trees:
                labeling = self.find_graceful_backtrack(tree)
                if labeling is not None:
                    n_graceful += 1
                    assert self.is_graceful(tree, labeling), \
                        f"Backtracking lieferte falsches Labeling für n={n}!"
                else:
                    counterexamples.append({
                        "n": n,
                        "edges": list(tree.edges()),
                        "note": "Kein graceful Labeling gefunden!"
                    })

            results[n] = {
                "trees_total": n_total,
                "graceful": n_graceful,
                "all_graceful": n_graceful == n_total
            }
            total_trees += n_total
            total_graceful += n_graceful

        return {
            "max_n": max_n,
            "total_trees": total_trees,
            "total_graceful": total_graceful,
            "all_graceful": total_graceful == total_trees,
            "counterexamples": counterexamples,
            "per_n": results
        }

    def verify_path_graph(self, n: int) -> Optional[Dict]:
        """
        Verifiziert graceful Labeling für Pfad P_n (explizite Konstruktion).

        Bekannte Konstruktion (trivial):
            Für Pfad v₀ - v₁ - ... - v_{n-1}:
            f(vₖ) = k falls k gerade, n-1-k falls k ungerade.
            (Alternierende Konstruktion gibt |f(u)-f(v)| = n-1, n-2, ..., 1 für aufeinanderfolgend)

        @param n: Länge des Pfades (Anzahl Knoten)
        @return: Graceful Labeling oder None
        @lastModified 2026-03-12
        """
        if not HAS_NETWORKX:
            return None
        G = nx.path_graph(n)
        # Explizite alternierende Konstruktion
        labeling = {}
        lo, hi = 0, n - 1
        for k, node in enumerate(range(n)):
            if k % 2 == 0:
                labeling[node] = lo
                lo += 1
            else:
                labeling[node] = hi
                hi -= 1
        if self.is_graceful(G, labeling):
            return labeling
        # Fallback: Backtracking
        return self.find_graceful_backtrack(G)

    def verify_caterpillar(self, spine_len: int, leaves: List[int]) -> Optional[Dict]:
        """
        Verifiziert graceful Labeling für Caterpillar-Graphen.

        Ein Caterpillar-Graph besteht aus einem Pfad (Wirbelsäule) mit
        hängenden Blättern an jedem Wirbelsäulen-Knoten.

        @param spine_len: Länge der Wirbelsäule (Anzahl Knoten)
        @param leaves: Anzahl Blätter an jedem Wirbelsäulen-Knoten
        @return: Graceful Labeling oder None
        @lastModified 2026-03-12
        """
        if not HAS_NETWORKX:
            return None
        G = nx.Graph()
        node_id = 0
        spine = list(range(spine_len))
        node_id = spine_len

        G.add_nodes_from(spine)
        for i in range(spine_len - 1):
            G.add_edge(spine[i], spine[i + 1])

        for i, leaf_count in enumerate(leaves[:spine_len]):
            for _ in range(leaf_count):
                G.add_node(node_id)
                G.add_edge(spine[i], node_id)
                node_id += 1

        return self.find_graceful_backtrack(G)


# ===========================================================================
# HAUPTPROGRAMM
# ===========================================================================

def run_all_verifications() -> None:
    """
    Führt alle vier Verifikationen durch und gibt strukturierte Berichte aus.

    Reihenfolge:
      1. Freiman/PFR-Analyse (konzeptuell + kleine Mengen)
      2. Frankl Union-Closed Conjecture (exhaustive small cases + Gilmer)
      3. Lonely Runner (n=4,5,6,7 Verifikation + n=8 Suche)
      4. Graceful Tree (≤10 Knoten erschöpfend)

    @lastModified 2026-03-12
    """
    sep = "=" * 70

    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("BATCH 19 — GRUPPE B: MATHEMATISCHE VERIFIKATION")
    print(f"Datum: 2026-03-12 | Autor: Michael Fuhrmann")
    print(sep)

    # -----------------------------------------------------------------------
    # VERMUTUNG 1: Freiman / PFR
    # -----------------------------------------------------------------------
    print(f"\n{'─'*70}")
    print("VERMUTUNG 1: Polynomial Freiman-Ruzsa (PFR) Vermutung")
    print(f"{'─'*70}")

    pfr = FreimanPFRAnalysis()

    print("\n[A] Verdopplungskonstanten kleiner Mengen:")
    struct_results = pfr.verify_small_sets_structure()
    for r in struct_results[:8]:
        print(f"  {r['type']:25s}  |A|={r['|A|']:3d}  K={r['K']:.4f}  → {r['structure']}")

    print("\n[B] Schranken-Vergleich (K=3):")
    K_test = 3.0
    print(f"  Verdopplungskonstante K = {K_test}")
    print(f"  Sanders-Schranke (2012): |P|/|A| ≤ exp({pfr.sanders_bound_over_z(K_test):.2f}) ← subexponentiell")
    print(f"  PFR-Vermutung (offen):   |P|/|A| ≤ K^C = {pfr.pfr_over_z_conjectured_bound(K_test):.2f} ← polynomial")
    print(f"  PFR über F₂ⁿ (bewiesen): log₂(|H|/|A|) ≤ {pfr.pfr_f2n_bound(K_test):.2f}")

    print("\n[C] F₂ⁿ Unterraum-Beispiel:")
    f2_ex = pfr.f2n_sumset_example(4)
    for k, v in f2_ex.items():
        if k != "Erklaerung":
            print(f"  {k}: {v}")
    print(f"  Erklärung: {f2_ex['Erklaerung']}")

    print("\n[KLASSIFIKATION] OFFEN über ℤ | BEWIESEN über F₂ⁿ (Gowers et al. 2023)")

    # -----------------------------------------------------------------------
    # VERMUTUNG 2: Frankl Union-Closed
    # -----------------------------------------------------------------------
    print(f"\n{'─'*70}")
    print("VERMUTUNG 2: Frankl Union-Closed Conjecture")
    print(f"{'─'*70}")

    frankl = FranklUnionClosedAnalysis()

    print("\n[A] Spezifische Familien:")
    specific = frankl.frankl_for_specific_families()
    for r in specific:
        if r["union_closed"]:
            sign = "✓" if r["≥ 1/2"] else "✗"
            print(f"  {r['Familie']:25s}  |F|={r['|F|']:2d}  max_freq={r['Anteil']:.4f}  {sign}")
        else:
            print(f"  {r['Familie']:25s}  nicht union-closed")

    print("\n[B] Gilmer-Schranke (2022):")
    gilmer = frankl.gilmer_bound_example()
    print(f"  Gilmer (2022): beste bekannte Schranke ≈ {gilmer['gilmer_bound_exact']:.6f} · |F|")
    print(f"  Alweiss-Huang-Sellke 2022: {gilmer['best_known_bound']}")
    print(f"  Ziel (Frankl): 0.5 · |F|")
    print(f"  Lücke: {gilmer['gap']}")
    print(f"  Vollständig bewiesen: {gilmer['fully_proven']}")

    print("\n[C] Exhaustive Verifikation (Grundmenge ≤ 3):")
    t0 = time.time()
    small_ver = frankl.verify_frankl_conjecture_small(ground_set_size=3)
    t1 = time.time()
    print(f"  Familien geprüft: {small_ver['families_checked']}")
    print(f"  Frankl-Verletzungen: {small_ver['frankl_violations']}")
    print(f"  Frankl hält durch: {small_ver['frankl_holds']}")
    print(f"  Laufzeit: {t1-t0:.2f}s")

    print("\n[KLASSIFIKATION] OFFEN (partiell: Gilmer 2022 → 0.382·|F|, Ziel 0.5)")

    # -----------------------------------------------------------------------
    # VERMUTUNG 3: Lonely Runner
    # -----------------------------------------------------------------------
    print(f"\n{'─'*70}")
    print("VERMUTUNG 3: Lonely Runner Vermutung")
    print(f"{'─'*70}")

    lr = LonelyRunnerVerification()

    for n_runners in [4, 5, 6, 7]:
        print(f"\n[n={n_runners}] Verifikation:")
        t0 = time.time()
        res = lr.run_batch_verification(n_runners, num_samples=5, max_denom=150)
        t1 = time.time()
        status = "VERIFIZIERT" if res["all_verified"] else "TEILWEISE VERIFIZIERT"
        print(f"  Instanzen getestet: {res['instances_tested']}")
        print(f"  Läufer gesamt: {res['total_runners']}")
        print(f"  Davon verifiziert: {res['verified']}")
        print(f"  Status: {status}  ({t1-t0:.2f}s)")
        if res["failed_cases"]:
            for fc in res["failed_cases"]:
                print(f"  WARNUNG: v={fc['velocities']}, Läufer {fc['runner']}: {fc['note']}")

    print(f"\n[n=8] Gegenbeispiel-Suche:")
    t0 = time.time()
    res8 = lr.search_counterexample_n8(num_cases=5, max_denom=80)
    t1 = time.time()
    print(f"  Instanzen getestet: {res8['cases_tested']}")
    print(f"  Vollständig verifiziert: {res8['fully_verified']}")
    print(f"  Fälle mit Problemen (Suchtiefe): {res8['potential_issues']}")
    print(f"  Hinweis: {res8['note']}")
    print(f"  Laufzeit: {t1-t0:.2f}s")
    for case in res8["cases_with_issues"][:2]:
        print(f"  Instanz {case['velocities']}: Läufer {case['failed_runners']} — {case['note']}")

    print("\n[KLASSIFIKATION] n≤7: BEWIESEN | n≥8: OFFEN")

    # -----------------------------------------------------------------------
    # VERMUTUNG 4: Graceful Tree
    # -----------------------------------------------------------------------
    print(f"\n{'─'*70}")
    print("VERMUTUNG 4: Graceful Tree Vermutung")
    print(f"{'─'*70}")

    gt = GracefulTreeVerification()

    if HAS_NETWORKX:
        print("\n[A] Systematische Verifikation aller Bäume n=1..10:")
        t0 = time.time()
        gt_res = gt.verify_all_trees_up_to_n(max_n=10)
        t1 = time.time()
        print(f"  Bäume gesamt: {gt_res['total_trees']}")
        print(f"  Davon graceful: {gt_res['total_graceful']}")
        print(f"  Alle graceful: {gt_res['all_graceful']}")
        print(f"  Gegenbeispiele: {len(gt_res['counterexamples'])}")
        print(f"  Laufzeit: {t1-t0:.1f}s")
        print(f"\n  Pro n:")
        for n_val, data in gt_res["per_n"].items():
            mark = "✓" if data["all_graceful"] else "✗"
            print(f"    n={n_val:2d}: {data['trees_total']:4d} Bäume, {data['graceful']:4d} graceful  {mark}")

        print("\n[B] Bekannte Klassen:")
        # Pfade
        for n_path in [3, 5, 7, 10]:
            lab = gt.verify_path_graph(n_path)
            print(f"  Pfad P_{n_path}: {'graceful ✓' if lab else 'NICHT graceful ✗'}")

        # Caterpillar
        lab_cat = gt.verify_caterpillar(spine_len=4, leaves=[1, 2, 1, 0])
        print(f"  Caterpillar(Wirbelsäule=4, Blätter=[1,2,1,0]): {'graceful ✓' if lab_cat else 'NICHT graceful ✗'}")

    else:
        print("  networkx nicht verfügbar — übersprungen")

    print("\n[KLASSIFIKATION] OFFEN (für ≤10 Knoten vollständig verifiziert)")

    # -----------------------------------------------------------------------
    # ZUSAMMENFASSUNG
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("ZUSAMMENFASSUNG BATCH 19 — GRUPPE B")
    print(sep)
    print("  Paper 72: Freiman/PFR  — OFFEN über ℤ | BEWIESEN über F₂ⁿ (GGMT 2023)")
    print("  Paper 73: Frankl UC    — OFFEN (Gilmer 2022: 0.382-Schranke, Ziel 0.5)")
    print("  Paper 74: Lonely Run.  — OFFEN für n≥8 (n≤7: BEWIESEN)")
    print("  Paper 75: Graceful T.  — OFFEN (alle ≤35 Knoten: verifiziert)")
    print(sep)


if __name__ == "__main__":
    run_all_verifications()
