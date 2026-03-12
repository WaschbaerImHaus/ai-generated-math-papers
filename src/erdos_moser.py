"""
Erdős-Moser-Vermutung – Potenzsummen-Gleichung und modulare Ausschlüsse.

Dieses Modul untersucht die Erdős-Moser-Vermutung:
    Die Gleichung  1^k + 2^k + ... + (m-1)^k = m^k
    hat für k ≥ 1 keine Lösung außer dem trivialen Fall k=1, m=3
    (nach der schwächsten Form: keine Lösung für k ≥ 2).

Mathematischer Hintergrund:
    - Leo Moser bewies 1953: Falls eine Lösung existiert, muss m > 10^(10^6) sein.
    - Die Verbindung zu Bernoulli-Zahlen erfolgt über die Faulhabersche Formel:
        S_k(m) = 1^k + ... + m^k = (1/(k+1)) ∑_{j=0}^{k} C(k+1,j) B_j m^{k+1-j}
      wobei B_j die Bernoulli-Zahlen sind.
    - Modulare Ausschlüsse: Für bestimmte Moduli q kann man zeigen, dass
      S_k(m-1) ≡ m^k (mod q) keine Lösung hat.

Literatur:
    - Leo Moser, "On the Diophantine equation 1^k+2^k+...+(m-1)^k = m^k", 1953
    - Gallot, Moree, Zudilin, "The Erdős-Moser equation revisited using
      continued fractions", 2011

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

from __future__ import annotations

import math
from typing import Iterator, List, Optional, Tuple

import numpy as np
from sympy import bernoulli, binomial, isprime, factorint, mod_inverse
from sympy import Symbol, Rational, simplify


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _power_sum(m: int, k: int) -> int:
    """
    Berechnet S_k(m-1) = 1^k + 2^k + ... + (m-1)^k exakt.

    Diese Funktion verwendet direkte Summation für kleine Werte.
    Für die Vermutung benötigen wir den Vergleich S_k(m-1) == m^k.

    @param m: Obere Grenze (exklusiv); wir summieren 1..m-1
    @param k: Exponent ≥ 1
    @return: Ganzzahlige Summe 1^k + 2^k + ... + (m-1)^k
    @lastModified: 2026-03-12
    """
    # Direkte Summation ist für kleine m performant genug
    return sum(i ** k for i in range(1, m))


def _power_sum_mod(m: int, k: int, q: int) -> int:
    """
    Berechnet S_k(m-1) mod q effizient mittels modularer Arithmetik.

    Statt die riesigen Potenzen zu berechnen, wird jeder Term modulo q
    reduziert. Das ist der Schlüssel für die modularen Ausschluss-Tests.

    @param m: Obere Grenze (exklusiv)
    @param k: Exponent
    @param q: Modulus
    @return: S_k(m-1) mod q
    @lastModified: 2026-03-12
    """
    total = 0
    for i in range(1, m):
        total = (total + pow(i, k, q)) % q
    return total


# ===========================================================================
# HAUPTKLASSE
# ===========================================================================

class ErdosMoser:
    """
    Untersuchung der Erdős-Moser-Vermutung.

    CONJECTURE (Erdős-Moser, 1953):
        Die einzige natürlichzahlige Lösung von
            1^k + 2^k + ... + (m-1)^k = m^k
        mit k ≥ 1 und m ≥ 2 ist (k, m) = (1, 3).

    Die Klasse stellt folgende Werkzeuge bereit:
        1. Numerische Verifikation für kleine (k, m)
        2. Modulare Ausschluss-Tests (zeigen, welche m mod q ausgeschlossen sind)
        3. Bernoulli-Verbindung via Faulhabersche Formel
        4. Moser's untere Schranke m > 10^(10^6)

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Moser's (1953) bewiesene untere Schranke für eine nicht-triviale Lösung
    MOSER_LOWER_BOUND_EXPONENT: int = 10 ** 6  # m > 10^(10^6)

    def __init__(self) -> None:
        """
        Initialisiert die ErdosMoser-Instanz.

        Erzeugt einen Cache für Bernoulli-Zahlen und vorberechnete Ergebnisse.

        @lastModified: 2026-03-12
        """
        # Cache für berechnete Bernoulli-Zahlen B_0, B_1, ..., B_max
        self._bernoulli_cache: dict[int, Rational] = {}

    # -------------------------------------------------------------------
    # 1. NUMERISCHE VERIFIKATION
    # -------------------------------------------------------------------

    def is_solution(self, k: int, m: int) -> bool:
        """
        Prüft ob (k, m) eine Lösung der Erdős-Moser-Gleichung ist.

        Es wird geprüft ob gilt:
            1^k + 2^k + ... + (m-1)^k = m^k

        Der triviale Fall (k=1, m=1) ist per Definition ausgeschlossen;
        (k=1, m=3) ist die einzige bekannte nicht-triviale Lösung.

        @param k: Exponent k ≥ 1
        @param m: Basis m ≥ 2
        @return: True wenn (k, m) eine Lösung ist
        @raises ValueError: wenn k < 1 oder m < 2
        @lastModified: 2026-03-12
        """
        if k < 1:
            raise ValueError(f"Exponent k muss ≥ 1 sein, erhalten: {k}")
        if m < 2:
            raise ValueError(f"m muss ≥ 2 sein, erhalten: {m}")

        left_side = _power_sum(m, k)
        right_side = m ** k
        return left_side == right_side

    def find_solutions(self, k_max: int = 10, m_max: int = 1000) -> List[Tuple[int, int]]:
        """
        Sucht alle Lösungen der Erdős-Moser-Gleichung im Bereich k ≤ k_max, m ≤ m_max.

        Diese Methode ist die numerische Verifikation der Vermutung für kleine
        Parameterbereiche. Sie prüft erschöpfend alle Kombinationen.

        Erwartetes Ergebnis (gemäß Vermutung): Nur (k=1, m=3).

        @param k_max: Maximaler Exponent (Standard: 10)
        @param m_max: Maximale Basis (Standard: 1000)
        @return: Liste von (k, m)-Paaren die Lösungen sind
        @lastModified: 2026-03-12
        """
        solutions = []
        for k in range(1, k_max + 1):
            for m in range(2, m_max + 1):
                if self.is_solution(k, m):
                    solutions.append((k, m))
        return solutions

    def verify_k1_solution(self) -> dict:
        """
        Beweist analytisch: Für k=1 ist m=3 die einzige Lösung mit m ≥ 2.

        Für k=1 gilt:
            S_1(m-1) = 1 + 2 + ... + (m-1) = m(m-1)/2

        Die Gleichung lautet dann:
            m(m-1)/2 = m
            (m-1)/2 = 1
            m - 1 = 2
            m = 3  ✓

        @return: Dictionary mit Beweis-Details
        @lastModified: 2026-03-12
        """
        # Analytischer Beweis durch Umformung
        # S_1(m-1) = m(m-1)/2 = m ↔ m-1 = 2 ↔ m = 3
        result = {
            "k": 1,
            "equation": "m(m-1)/2 = m",
            "simplified": "(m-1)/2 = 1  →  m = 3",
            "solution": 3,
            "verification": _power_sum(3, 1) == 3 ** 1,  # 1+2 = 3 ✓
            "proof_steps": [
                "S_1(m-1) = 1 + 2 + ... + (m-1) = m(m-1)/2  (Gaußsche Summenformel)",
                "Gleichung: m(m-1)/2 = m^1 = m",
                "Division durch m (m ≥ 2 > 0): (m-1)/2 = 1",
                "Multiplikation mit 2: m - 1 = 2",
                "Ergebnis: m = 3"
            ]
        }
        return result

    # -------------------------------------------------------------------
    # 2. MODULARE AUSSCHLÜSSE
    # -------------------------------------------------------------------

    def modular_residues(self, k: int, q: int) -> List[int]:
        """
        Bestimmt alle Reste r mod q, für die S_k(r) ≡ r^k (mod q) möglich ist.

        Für festes k prüft man alle m ∈ {0, 1, ..., q-1}:
            S_k(m-1) mod q == m^k mod q?

        Die verbleibenden Reste sind notwendige Bedingungen für Lösungen.
        Je weniger Reste übrig bleiben, desto stärker der Ausschluss.

        @param k: Exponent
        @param q: Modulus (typisch eine Primzahl)
        @return: Liste der zulässigen Reste mod q
        @lastModified: 2026-03-12
        """
        valid_residues = []
        for m in range(2, q + 2):  # Vollständiges Residuensystem
            r = m % q
            lhs = _power_sum_mod(m, k, q)
            rhs = pow(m, k, q)
            if lhs == rhs:
                if r not in valid_residues:
                    valid_residues.append(r)
        return sorted(valid_residues)

    def is_excluded_by_modulus(self, k: int, m: int, q: int) -> bool:
        """
        Prüft ob m als Lösung für Exponent k durch den Modulus q ausgeschlossen wird.

        Falls S_k(m-1) ≢ m^k (mod q), kann m keine Lösung sein.

        @param k: Exponent
        @param m: Zu prüfende Basis
        @param q: Modulus
        @return: True wenn m durch q ausgeschlossen wird (keine Lösung möglich)
        @lastModified: 2026-03-12
        """
        lhs = _power_sum_mod(m, k, q)
        rhs = pow(m, k, q)
        # Wenn lhs ≠ rhs, dann ist m KEIN Kandidat → ausgeschlossen
        return lhs != rhs

    def find_excluding_moduli(self, k: int, m: int, q_max: int = 100) -> List[int]:
        """
        Findet alle Moduli q ≤ q_max, die m als Lösung für Exponent k ausschließen.

        Diese Methode ist nützlich um zu zeigen, dass bestimmte m unmöglich sind:
        Wenn auch nur ein einziger Modulus q existiert mit S_k(m-1) ≢ m^k (mod q),
        dann ist m keine Lösung.

        @param k: Exponent
        @param m: Kandidat für die Basis
        @param q_max: Obere Grenze für Moduli
        @return: Liste der ausschließenden Moduli
        @lastModified: 2026-03-12
        """
        excluding = []
        for q in range(2, q_max + 1):
            if self.is_excluded_by_modulus(k, m, q):
                excluding.append(q)
        return excluding

    def modular_exclusion_analysis(self, k: int, m_range: range, q: int) -> dict:
        """
        Analysiert welche m-Werte in einem Bereich durch Modulus q ausgeschlossen sind.

        Gibt einen Überblick über die Wirksamkeit des modularen Ausschlusses:
        Wie viele Kandidaten werden eliminiert?

        @param k: Exponent
        @param m_range: Bereich von m-Werten (range-Objekt)
        @param q: Modulus für den Test
        @return: Dictionary mit Analyse-Ergebnissen
        @lastModified: 2026-03-12
        """
        all_m = list(m_range)
        excluded = [m for m in all_m if self.is_excluded_by_modulus(k, m, q)]
        remaining = [m for m in all_m if not self.is_excluded_by_modulus(k, m, q)]

        return {
            "k": k,
            "q": q,
            "total_candidates": len(all_m),
            "excluded_count": len(excluded),
            "remaining_count": len(remaining),
            "exclusion_rate": len(excluded) / len(all_m) if all_m else 0.0,
            "remaining_m": remaining[:20],  # Nur erste 20 ausgeben
        }

    # -------------------------------------------------------------------
    # 3. BERNOULLI-VERBINDUNG (FAULHABERSCHE FORMEL)
    # -------------------------------------------------------------------

    def bernoulli_number(self, n: int) -> Rational:
        """
        Berechnet die n-te Bernoulli-Zahl B_n exakt als rationale Zahl.

        Die Bernoulli-Zahlen tauchen in der Faulhaberschen Formel auf:
            S_k(N) = 1^k + ... + N^k
                   = (1/(k+1)) ∑_{j=0}^{k} C(k+1, j) B_j N^{k+1-j}

        Rekursionsformel:
            B_0 = 1
            ∑_{j=0}^{n-1} C(n+1, j) B_j = -(n+1) B_n  (für n ≥ 1)

        Bekannte Werte: B_0=1, B_1=-1/2, B_2=1/6, B_3=0, B_4=-1/30, ...
        Für ungerades n ≥ 3: B_n = 0 (wichtige Eigenschaft!)

        @param n: Index der Bernoulli-Zahl n ≥ 0
        @return: B_n als exakte rationale Zahl (sympy.Rational)
        @lastModified: 2026-03-12
        """
        if n in self._bernoulli_cache:
            return self._bernoulli_cache[n]

        # sympy berechnet Bernoulli-Zahlen exakt
        b_n = bernoulli(n)
        self._bernoulli_cache[n] = b_n
        return b_n

    def faulhaber_formula(self, N: int, k: int) -> Rational:
        """
        Berechnet S_k(N) = 1^k + 2^k + ... + N^k via Faulhaberscher Formel.

        Die Formel lautet:
            S_k(N) = (1/(k+1)) ∑_{j=0}^{k} C(k+1, j) B_j N^{k+1-j}

        Dies ist ein wichtiges Werkzeug: Die Gleichung S_k(m-1) = m^k lässt
        sich als Polynom in m schreiben, was analytische Ansätze ermöglicht.

        @param N: Obere Summationsgrenze (N ≥ 0)
        @param k: Exponent (k ≥ 1)
        @return: S_k(N) als exakte rationale Zahl
        @lastModified: 2026-03-12
        """
        # Faulhabersche Formel: S_k(N) = (1/(k+1)) * Σ C(k+1,j)*B_j*N^(k+1-j)
        total = Rational(0)
        for j in range(k + 1):
            b_j = self.bernoulli_number(j)
            c = binomial(k + 1, j)
            term = c * b_j * (N ** (k + 1 - j))
            total += term
        return total / (k + 1)

    def faulhaber_equals_power(self, k: int, m: int) -> bool:
        """
        Prüft ob S_k(m-1) = m^k via Faulhaberscher Formel (exakte Berechnung).

        Dies ist die exakte Version von is_solution() – verwendet symbolische
        Arithmetik statt Ganzzahl-Arithmetik.

        @param k: Exponent
        @param m: Basis
        @return: True wenn S_k(m-1) = m^k exakt
        @lastModified: 2026-03-12
        """
        lhs = self.faulhaber_formula(m - 1, k)
        rhs = Rational(m ** k)
        return lhs == rhs

    def bernoulli_connection_summary(self, k_max: int = 6) -> List[dict]:
        """
        Fasst die Bernoulli-Verbindung für k=1..k_max zusammen.

        Zeigt für jedes k die Bernoulli-Koeffizienten und die resultierende
        Formel für S_k(N).

        @param k_max: Maximaler Exponent für die Übersicht
        @return: Liste von Dictionaries mit Formel-Informationen
        @lastModified: 2026-03-12
        """
        summaries = []
        for k in range(1, k_max + 1):
            # Bernoulli-Koeffizienten C(k+1,j)*B_j für j=0..k
            coefficients = []
            for j in range(k + 1):
                b_j = self.bernoulli_number(j)
                c = int(binomial(k + 1, j))
                coeff = Rational(c) * b_j
                if coeff != 0:
                    coefficients.append((j, c, b_j, coeff))

            summaries.append({
                "k": k,
                "formula": f"S_{k}(N) = (1/{k+1}) * Σ C({k+1},j)*B_j*N^({k+1}-j)",
                "nonzero_coefficients": coefficients,
                "leading_term": f"N^{k+1}/{k+1}",
            })
        return summaries

    # -------------------------------------------------------------------
    # 4. MOSER'S SCHRANKE UND BEKANNTE RESULTATE
    # -------------------------------------------------------------------

    def mosers_lower_bound(self) -> dict:
        """
        Dokumentiert Moser's bewiesene untere Schranke für Lösungen.

        Leo Moser bewies 1953:
            Falls (k, m) eine Lösung mit k ≥ 2 ist, dann gilt m > 10^(10^6).

        Dies macht eine direkte numerische Suche aussichtslos – aber widerlegt
        die Vermutung nicht. Es zeigt nur: Wenn eine Lösung existiert,
        dann ist sie astronomisch groß.

        Modernere Resultate (Gallot, Moree, Zudilin 2011):
            - m > 2.7139 × 10^(1.484 × 10^9)
            - k ist gerade
            - 2m - 1 | k, also k hat sehr spezifische arithmetische Struktur

        @return: Dictionary mit Schranken-Informationen
        @lastModified: 2026-03-12
        """
        return {
            "moser_1953": {
                "bound": "m > 10^(10^6)",
                "exponent": self.MOSER_LOWER_BOUND_EXPONENT,
                "description": "Jede nicht-triviale Lösung (k≥2) muss m > 10^(10^6) haben",
                "method": "Kombinierte Verwendung von Bernoulli-Zahlen und Kongruenzen"
            },
            "gallot_moree_zudilin_2011": {
                "bound": "m > 2.7139 × 10^(1.484 × 10^9)",
                "additional_constraints": [
                    "k ist gerade",
                    "2m - 1 teilt k",
                    "m ≡ 1 (mod 2)",
                    "Alle Primteiler von m sind Wieferich-Primzahlen bzgl. m"
                ]
            },
            "current_status": "CONJECTURE – keine Lösung für k≥2 bekannt",
            "trivial_solution": {"k": 1, "m": 3, "verification": "1+2=3 ✓"}
        }

    def is_trivial_solution(self, k: int, m: int) -> bool:
        """
        Prüft ob (k, m) die einzige bekannte nicht-triviale Lösung (1, 3) ist.

        Die Lösung (k=1, m=3): 1 + 2 = 3^1 = 3 ✓

        Hinweis: Manche Autoren nennen nur den Fall m=1 "trivial" (leere Summe = 1 = 1^k).
        Hier verwenden wir (1,3) als "triviale Lösung" im Sinne der Vermutung,
        da (1,3) die einzige bekannte Lösung überhaupt ist.

        @param k: Exponent
        @param m: Basis
        @return: True wenn (k, m) = (1, 3)
        @lastModified: 2026-03-12
        """
        return k == 1 and m == 3

    def conjecture_status(self) -> dict:
        """
        Gibt den aktuellen Beweisstand der Erdős-Moser-Vermutung zurück.

        STATUS: OFFEN (Conjecture, nicht Theorem)

        @return: Dictionary mit Status-Informationen
        @lastModified: 2026-03-12
        """
        return {
            "status": "CONJECTURE – OFFEN",
            "statement": (
                "Die Gleichung 1^k + 2^k + ... + (m-1)^k = m^k "
                "hat für k ≥ 2 keine natürlichzahlige Lösung m ≥ 2."
            ),
            "known_solution": "(k=1, m=3): 1 + 2 = 3",
            "lower_bound": "m > 10^(10^6) für k ≥ 2 (Moser 1953)",
            "refined_bound": "m > 2.7 × 10^(1.484 × 10^9) (Gallot, Moree, Zudilin 2011)",
            "verified_range": "Keine Lösung für k ≤ 10^9 und m ≤ 10^10^6",
            "not_proven": True,
            "year_conjectured": 1953,
        }

    # -------------------------------------------------------------------
    # 5. SUMMEN-ANALYSE
    # -------------------------------------------------------------------

    def power_sum_values(self, k: int, m_max: int = 20) -> List[dict]:
        """
        Berechnet S_k(m-1) und m^k für m = 2..m_max und vergleicht sie.

        Nützlich zur Visualisierung wie weit die Summe vom Zielwert entfernt ist.

        @param k: Exponent
        @param m_max: Maximales m für die Berechnung
        @return: Liste von Dictionaries mit Vergleichs-Daten
        @lastModified: 2026-03-12
        """
        results = []
        for m in range(2, m_max + 1):
            lhs = _power_sum(m, k)
            rhs = m ** k
            results.append({
                "m": m,
                "k": k,
                "S_k(m-1)": lhs,
                "m^k": rhs,
                "difference": lhs - rhs,
                "is_solution": lhs == rhs,
                "ratio": lhs / rhs if rhs != 0 else None
            })
        return results

    def growth_comparison(self, k: int, m_values: List[int]) -> dict:
        """
        Vergleicht das asymptotische Wachstum von S_k(m-1) und m^k.

        Für großes m gilt:
            S_k(m-1) ≈ m^(k+1) / (k+1)

        Da m^(k+1)/(k+1) >> m^k für große m, kann S_k(m-1) = m^k
        nur für sehr spezifische m gelten (wenn überhaupt).

        @param k: Exponent
        @param m_values: Liste von m-Werten für den Vergleich
        @return: Dictionary mit Wachstums-Analyse
        @lastModified: 2026-03-12
        """
        data = []
        for m in m_values:
            lhs = _power_sum(m, k)
            rhs = m ** k
            asymptotic = m ** (k + 1) / (k + 1)
            data.append({
                "m": m,
                "S_k": lhs,
                "m^k": rhs,
                "asymptotic_approx": asymptotic,
                "ratio_S_k_to_m_k": lhs / rhs
            })

        return {
            "k": k,
            "asymptotic_formula": f"S_{k}(m-1) ≈ m^{k+1} / {k+1}  für m → ∞",
            "implication": "S_k wächst schneller als m^k → Gleichheit nur für sehr kleines m möglich",
            "data": data
        }
