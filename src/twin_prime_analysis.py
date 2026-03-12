"""
@file twin_prime_analysis.py
@brief Analyse von Zwillingsprimzahlen: Statistik, Hardy-Littlewood, Zhang/Maynard.
@description
    Dieses Modul untersucht Zwillingsprimzahlpaare (p, p+2) umfassend:
    Zählung, Dichte, Lückenstatistik und moderne Durchbrüche.

    **Zwillingsprimzahl-Vermutung** (Conjecture, offen):
        Es gibt unendlich viele Primzahlpaare (p, p+2).
        Bekannt: (3,5), (5,7), (11,13), (17,19), (29,31), ...

    **Hardy-Littlewood Vermutung B** (1923, Conjecture):
        π₂(x) ~ 2·C₂·x / (ln x)²
        mit C₂ = ∏_{p>2} p(p−2)/(p−1)² ≈ 0.6601618158...

    **Brun-Konstante** B₂ (Viggo Brun, 1919 — BEWEIS):
        B₂ = Σ_{(p,p+2) twin prime} (1/p + 1/(p+2)) konvergiert ≈ 1.9021605831...

    **Moderne Durchbrüche**:
    - Zhang (2013): Es gibt unendlich viele Primpaare mit Abstand < 70.000.000
    - Maynard (2013): Abstand < 600 (GPY-Methode verbessert)
    - Polymath8b (2014): Abstand < 246
    - Vermutung: Abstand = 2 (Zwillingsprimzahl-Vermutung) — OFFEN

    **GPY-Sieb-Technik** (Goldston-Pintz-Yıldırım 2005):
        Zeigt bedingte Häufungen von Primzahlen, aber nicht Abstand 2.
        Basis: Siebbindungsmethoden + Mollifier-Technik.

    **Parity Problem** (Selberg):
        Siebmethoden können Primzahlpaare (p, p+2) nicht direkt zählen,
        weil die Paritätsbeschränkung den Sieb "blind" für Abstand 2 macht.

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np


# ===========================================================================
# HILFSFUNKTIONEN: Siebe
# ===========================================================================

def _sieve_of_eratosthenes(limit: int) -> np.ndarray:
    """
    Erathostenes-Sieb: Berechnet alle Primzahlen bis limit.

    Zeitkomplexität: O(n log log n)
    Speicher: O(n) Bits

    @param limit: Obere Grenze (inklusiv)
    @return: Boolean-Array, is_prime[i] = True wenn i prim
    @lastModified: 2026-03-12
    """
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if is_prime[i]:
            is_prime[i * i::i] = False
    return is_prime


# ===========================================================================
# HAUPTKLASSE: TwinPrimeAnalysis
# ===========================================================================

class TwinPrimeAnalysis:
    """
    Analyse von Zwillingsprimzahlpaaren und verwandten Primzahllücken.

    Berechnet alle Zwillingsprimpaare bis zu einem Limit und analysiert
    deren Dichte, Verteilung und Vergleich mit Hardy-Littlewood-Heuristik.

    **Hardy-Littlewood Konstante C₂** (Zwillingsprimzahl-Konstante):
        C₂ = ∏_{p≥3, p prim} p(p−2)/(p−1)² ≈ 0.6601618...

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    # Hardy-Littlewood Zwillingsprimzahl-Konstante C₂
    # Berechnet als Grenzwert des Eulerprodukts ∏_{p≥3} p(p−2)/(p−1)²
    C2_HARDY_LITTLEWOOD = 0.6601618158468695739278121

    def __init__(self, limit: int = 100_000):
        """
        Initialisiert die Analyse bis zum gegebenen Limit.

        @param limit: Obere Schranke für Primzahlsuche
        @lastModified: 2026-03-12
        """
        self.limit = limit
        self._is_prime: Optional[np.ndarray] = None
        self._twin_pairs: Optional[List[Tuple[int, int]]] = None

    def _ensure_sieve(self) -> None:
        """
        Berechnet den Sieb wenn noch nicht vorhanden.
        @lastModified: 2026-03-12
        """
        if self._is_prime is None:
            self._is_prime = _sieve_of_eratosthenes(self.limit)

    def get_twin_pairs(self) -> List[Tuple[int, int]]:
        """
        Gibt alle Zwillingsprimpaare (p, p+2) bis self.limit zurück.

        Ein Zwillingsprimpaar (p, p+2) besteht aus zwei Primzahlen
        deren Differenz genau 2 beträgt.

        @return: Liste aller Zwillingsprimpaare (p, p+2) mit p ≤ limit−2
        @lastModified: 2026-03-12
        """
        if self._twin_pairs is not None:
            return self._twin_pairs
        self._ensure_sieve()
        ip = self._is_prime
        # Alle p mit ip[p] und ip[p+2] True
        primes = np.where(ip)[0]
        twins = []
        for p in primes:
            if p + 2 <= self.limit and ip[p + 2]:
                twins.append((int(p), int(p + 2)))
        self._twin_pairs = twins
        return twins

    def count_twin_primes(self, x: Optional[int] = None) -> int:
        """
        Gibt die Anzahl der Zwillingsprimpaare bis x (oder self.limit) zurück.

        π₂(x) = #{(p, p+2) : p ≤ x, p und p+2 prim}

        @param x: Obere Schranke (Standard: self.limit)
        @return: Anzahl der Zwillingsprimpaare
        @lastModified: 2026-03-12
        """
        pairs = self.get_twin_pairs()
        if x is None:
            return len(pairs)
        return sum(1 for p, _ in pairs if p <= x)

    def hardy_littlewood_prediction(self, x: int) -> float:
        """
        Hardy-Littlewood-Vorhersage für die Anzahl der Zwillingsprimpaare bis x.

        **Vermutung** (Hardy-Littlewood 1923, Conjecture B — UNBEWIESEN):
            π₂(x) ~ 2·C₂·x / (ln x)²

        @param x: Obere Schranke x ≥ 2
        @return: Vorhergesagte Anzahl der Zwillingsprimpaare
        @lastModified: 2026-03-12
        """
        if x < 4:
            return 0.0
        return 2.0 * self.C2_HARDY_LITTLEWOOD * x / (math.log(x) ** 2)

    def compute_C2_numerically(self, prime_limit: int = 10_000) -> float:
        """
        Berechnet die Hardy-Littlewood-Konstante C₂ numerisch als Eulerprodukt.

        C₂ = ∏_{p≥3, p prim} p(p−2)/(p−1)²

        Das Produkt konvergiert sehr langsam. Für bessere Genauigkeit werden
        mehr Primzahlen benötigt.

        @param prime_limit: Obere Schranke für das Eulerprodukt
        @return: Numerische Näherung von C₂
        @lastModified: 2026-03-12
        """
        is_prime = _sieve_of_eratosthenes(prime_limit)
        primes = np.where(is_prime)[0]
        # Starte mit p=3 (p=2 ist separat behandelt, Faktor = 1 da 2(2-2)/(2-1)² = 0)
        # Korrektes Eulerprodukt: Faktor für p=2 ist 2/(1)² nicht enthalten
        product = 1.0
        for p in primes:
            if p <= 2:
                continue
            factor = p * (p - 2) / ((p - 1) ** 2)
            product *= factor
        return product

    def prime_gaps_of_size(self, gap: int) -> List[int]:
        """
        Gibt alle Primzahlen p zurück, bei denen p_{next} - p = gap.

        Für gap=2: Zwillingsprimzahlen.
        Für gap=4: Cousin-Primzahlen.
        Für gap=6: Sexy-Primzahlen.

        @param gap: Primzahllückengröße (positiv, gerade für p>2)
        @return: Liste der kleineren Primzahlen der Paare
        @lastModified: 2026-03-12
        """
        self._ensure_sieve()
        ip = self._is_prime
        primes = list(np.where(ip)[0])
        result = []
        for i in range(len(primes) - 1):
            if primes[i + 1] - primes[i] == gap:
                result.append(int(primes[i]))
        return result

    def gap_statistics(self, max_gap: int = 20) -> Dict[int, int]:
        """
        Berechnet die Häufigkeit jeder Primzahllückengröße bis max_gap.

        @param max_gap: Maximale Lückengröße
        @return: Dictionary {gap: count}
        @lastModified: 2026-03-12
        """
        self._ensure_sieve()
        primes = list(np.where(self._is_prime)[0])
        counts: Dict[int, int] = {}
        for i in range(len(primes) - 1):
            g = primes[i + 1] - primes[i]
            if g <= max_gap:
                counts[g] = counts.get(g, 0) + 1
        return counts

    def cramer_model_prediction(self, x: int, gap: int = 2) -> float:
        """
        Cramér-Modell-Vorhersage für Primzahllücken der Größe gap.

        Cramér-Modell (1936, heuristisch):
            Jede Zahl n ist "unabhängig" prim mit Wahrscheinlichkeit 1/ln(n).
            → Erwartete Anzahl von Paaren (p, p+gap) bis x:
                E[pairs] ≈ x / (ln x)²  (für gap=2 mit C₂-Korrekturfaktor)

        @param x: Obere Schranke
        @param gap: Lückengröße
        @return: Cramér-Modell-Vorhersage
        @lastModified: 2026-03-12
        """
        if x < 4:
            return 0.0
        # Einfaches Cramér-Modell ohne Siebkorrekturen
        return x / (math.log(x) ** 2)

    def zhang_maynard_summary(self) -> Dict:
        """
        Zusammenfassung der modernen Durchbrüche bei begrenzten Primzahllücken.

        Zhang (2013), Maynard (2013), Polymath8b (2014) zeigten:
        Es gibt unendlich viele Primzahlpaare mit beschränktem Abstand.

        @return: Dictionary mit historischen Ergebnissen
        @lastModified: 2026-03-12
        """
        return {
            "Zwillingsprimzahl-Vermutung (offen)": {
                "bound": 2,
                "status": "Conjecture — unbewiesen",
                "year": "offen",
            },
            "Zhang (2013)": {
                "bound": 70_000_000,
                "status": "Beweis: lim inf (p_{n+1} - p_n) < 7×10^7",
                "year": 2013,
                "method": "GPY-Sieb + Siebbindung + Zhang-Verfeinerung",
            },
            "Maynard (2013)": {
                "bound": 600,
                "status": "Beweis: lim inf (p_{n+1} - p_n) < 600",
                "year": 2013,
                "method": "Gewichtetes GPY-Sieb (m-Tupel-Variante)",
            },
            "Polymath8b (2014)": {
                "bound": 246,
                "status": "Beweis: lim inf (p_{n+1} - p_n) ≤ 246",
                "year": 2014,
                "method": "Optimierte Maynard-Siebe + kollaborative Optimierung",
            },
        }

    def gpy_sieve_explanation(self) -> str:
        """
        Erklärt die GPY-Siebtechnik (Goldston-Pintz-Yıldırım).

        GPY (2005) zeigt: Für jedes ε > 0 gibt es Primzahlpaare mit Abstand < ε·ln p.
        Das bedeutet: Primzahlen sind "viel häufiger nahe beieinander" als erwartet.
        Aber der Schritt zu einem festen Abstand (wie 2) gelingt nicht.

        @return: Erklärungstext
        @lastModified: 2026-03-12
        """
        return (
            "GPY-Siebtechnik (Goldston-Pintz-Yıldırım, 2005):\n"
            "Grundidee: Betrachte gewichtete Summe S = Σ_{n~N} w(n) · Λ(n+h₁)·Λ(n+h₂)\n"
            "mit Mollifier w(n) = (Σ_{d|P(n)} λ_d)², P(n) = n(n+2).\n"
            "GPY zeigt: S > 0 für geeignete λ_d, also gibt es Paare (n+h₁, n+h₂)\n"
            "mit beiden Primzahlen und h₁−h₂ beliebig klein (aber > 2).\n\n"
            "Parity Problem (Selberg):\n"
            "Das Parity Problem verhindert, dass klassische Siebe direkt\n"
            "Paare mit Abstand 2 zählen. Das Sieb kann nicht zwischen\n"
            "p(p+2) und p·q·(p+2)·r unterscheiden, wenn beide die gleiche\n"
            "Siebstruktur haben (beide haben 0 oder 2 Primteiler ≤ √x).\n"
            "Zhang umgeht dies durch verschärfte Exponentialsummen-Schranken."
        )

    def selberg_parity_problem(self) -> str:
        """
        Erläutert das Selberg-Parity-Problem für Zwillingsprimzahlsiebe.

        @return: Erklärungstext
        @lastModified: 2026-03-12
        """
        return (
            "Selberg-Parity-Problem:\n"
            "Sei f(n) die Anzahl der Primfaktoren von n (mit Vielfachheit).\n"
            "Dann ist f(n) + f(n+2) mod 2 für Siebfunktionen nicht zu kontrollieren.\n"
            "Für Zwillingsprimpaare: f(n) = 1 und f(n+2) = 1 → Summe = 2 (gerade).\n"
            "Für n = p·q, n+2 = r: f(n) = 2, f(n+2) = 1 → Summe = 3 (ungerade).\n"
            "Klassische Siebe können diese beiden Fälle nicht trennen,\n"
            "weil die Siebfunktion Λ²(n) auf beiden den gleichen Wert hat.\n"
            "→ Man kann daher nicht Σ Λ(n)Λ(n+2) korrekt abschätzen."
        )


# ===========================================================================
# KLASSE: PrimeSieveGaps
# ===========================================================================

class PrimeSieveGaps:
    """
    Analyse von Primzahllücken: Cramér-Modell vs. tatsächliche Verteilung,
    sowie Erdős-Rankin-Konstruktion für große Lücken.

    **Cramér-Vermutung** (Harald Cramér, 1936, Conjecture — UNBEWIESEN):
        lim sup_{n→∞} (p_{n+1} − p_n) / (ln p_n)² = 1

    **Erdős-Rankin-Konstruktion** (Erdős 1935, Rankin 1938):
        Es gibt unendlich viele Primzahllücken der Größe
        g_n > c · (ln p_n · ln ln p_n · ln ln ln ln p_n) / (ln ln ln p_n)²

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    def __init__(self, limit: int = 1_000_000):
        """
        Initialisiert die Lückenanalyse bis limit.

        @param limit: Obere Schranke für Primzahlsuche
        @lastModified: 2026-03-12
        """
        self.limit = limit
        self._primes: Optional[List[int]] = None

    def _ensure_primes(self) -> None:
        """Berechnet Primzahlliste wenn noch nicht vorhanden."""
        if self._primes is None:
            is_prime = _sieve_of_eratosthenes(self.limit)
            self._primes = list(map(int, np.where(is_prime)[0]))

    def get_all_gaps(self) -> List[Tuple[int, int]]:
        """
        Gibt alle Primzahllücken (p_n, gap_n) zurück.

        @return: Liste von (p_n, p_{n+1} - p_n) für alle Primzahlen bis limit
        @lastModified: 2026-03-12
        """
        self._ensure_primes()
        primes = self._primes
        return [(primes[i], primes[i + 1] - primes[i]) for i in range(len(primes) - 1)]

    def gaps_by_size(self, size: int) -> List[int]:
        """
        Gibt alle Primzahlen zurück, nach denen eine Lücke der Größe 'size' folgt.

        @param size: Lückengröße
        @return: Liste der Primzahlen p mit p_{next} - p = size
        @lastModified: 2026-03-12
        """
        return [p for p, g in self.get_all_gaps() if g == size]

    def cramers_model_comparison(self) -> Dict:
        """
        Vergleicht tatsächliche Lückenverteilung mit Cramér-Modell-Vorhersage.

        Cramér-Modell: Primzahlen sind "lokal Poisson-verteilt" mit Rate 1/ln(p).
        → Lücken folgen Exponentialverteilung mit Parameter 1/ln(p).

        @return: Vergleichs-Dictionary
        @lastModified: 2026-03-12
        """
        self._ensure_primes()
        primes = self._primes
        if len(primes) < 2:
            return {}

        gaps = [primes[i + 1] - primes[i] for i in range(len(primes) - 1)]

        # Tatsächliche Statistik
        actual_mean = sum(gaps) / len(gaps)
        actual_max = max(gaps)

        # Cramér-Vorhersage: mittlere Lücke ≈ ln(N) am Ende des Intervalls
        cramer_mean = math.log(self.limit)

        # Standardabweichung
        variance = sum((g - actual_mean) ** 2 for g in gaps) / len(gaps)
        std_dev = math.sqrt(variance)

        return {
            "total_gaps": len(gaps),
            "actual_mean_gap": actual_mean,
            "cramer_predicted_mean": cramer_mean,
            "actual_max_gap": actual_max,
            "cramer_max_prediction": math.log(self.limit) ** 2,
            "std_dev": std_dev,
            "gaps_of_size_2": gaps.count(2),
            "gaps_of_size_4": gaps.count(4),
            "gaps_of_size_6": gaps.count(6),
        }

    def erdos_rankin_lower_bound(self, x: float) -> float:
        """
        Erdős-Rankin-Konstruktion: untere Schranke für maximale Primzahllücke bis x.

        Satz (Erdős 1935, Rankin 1938):
            Es gibt Primzahllücken der Größe ≥ c·ln(x)·ln ln(x)·ln ln ln ln(x) / (ln ln ln(x))²

        Dies ist eine BEWIESENE untere Schranke (kein Rekordwert, nur Existenzaussage).

        @param x: Obere Schranke x > ee^e (= e^{e^e} ≈ 5.3·10^6)
        @return: Untere Schranke für maximale Primzahllücke
        @lastModified: 2026-03-12
        """
        if x <= math.e:
            return 0.0
        ln_x = math.log(x)
        ln_ln_x = math.log(ln_x) if ln_x > 1 else 1.0
        ln_ln_ln_x = math.log(ln_ln_x) if ln_ln_x > 1 else 1.0
        ln_ln_ln_ln_x = math.log(ln_ln_ln_x) if ln_ln_ln_x > 1 else 1.0
        # c ≈ 1 (Rankin-Konstante, verbessert durch Maier 1981 auf c = e^γ)
        c = math.exp(0.5772156649)  # e^γ (γ = Euler-Mascheroni-Konstante)
        numerator = ln_x * ln_ln_x * ln_ln_ln_ln_x
        denominator = ln_ln_ln_x ** 2
        return c * numerator / denominator if denominator > 0 else 0.0

    def maximal_gap_record(self) -> Tuple[int, int, int]:
        """
        Findet die größte Primzahllücke bis self.limit.

        @return: Tupel (p, p_next, gap) der maximalen Lücke
        @lastModified: 2026-03-12
        """
        self._ensure_primes()
        primes = self._primes
        if len(primes) < 2:
            return (2, 3, 1)
        max_gap = 0
        max_p = primes[0]
        max_p_next = primes[1]
        for i in range(len(primes) - 1):
            g = primes[i + 1] - primes[i]
            if g > max_gap:
                max_gap = g
                max_p = primes[i]
                max_p_next = primes[i + 1]
        return (max_p, max_p_next, max_gap)
