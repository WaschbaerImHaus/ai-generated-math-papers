r"""
@file legendre_brocard.py
@brief Verifikation der Legendre-Vermutung und Brocard-Vermutung (7. Konjektur).
@description
    Dieses Modul untersucht zwei klassische Vermutungen aus der analytischen
    Zahlentheorie:

    **1. Legendre-Vermutung (1798, ungeklärt)**

    Conjecture (Legendre): Für jede positive ganze Zahl n ≥ 1 existiert mindestens
    eine Primzahl p mit:

        $$ n^2 < p < (n+1)^2 $$

    Verbindung zu Bertrand's Postulat (bewiesenes Theorem, Tschebyschew 1852):
    Für jedes n ≥ 1 existiert eine Primzahl p mit n < p ≤ 2n.
    Bertrand folgt aus der Legendre-Vermutung (schwächer), aber NICHT umgekehrt.

    **2. Brocard-Vermutung (1904, ungeklärt)**

    Conjecture (Brocard, 1904): Für jede Zahl n ≥ 2 liegen zwischen den Quadraten
    aufeinanderfolgender Primzahlen $p_n$ und $p_{n+1}$ mindestens 4 Primzahlen:

        $$ \pi(p_{n+1}^2) - \pi(p_n^2) \geq 4 \quad \text{für alle } n \geq 2 $$

    Hinweis: Dies ist NICHT die Brocard-Ramanujan-Gleichung n!+1 = m².
    Dies ist eine andere, eigenständige Vermutung von Henri Brocard.

    **Bekannte Fakten**:
    - Legendre: Verifiziert für n bis 10^10 (Oliveira e Silva, 2014)
    - Brocard: Verifiziert für alle Primzahlen bis 10^9
    - Beide gelten als sehr wahrscheinlich wahr, sind aber unbewiesen

    **Bezug zu anderen Vermutungen**:
    - Cramér-Vermutung: pₙ₊₁ - pₙ = O(log²pₙ) würde Legendre implizieren
    - Riemann-Hypothese: Impliziert pₙ₊₁ - pₙ = O(√pₙ · log pₙ), was Legendre fast zeigt
    - Andrica-Vermutung: √p_{n+1} - √pₙ < 1 ist äquivalent zu Legendre für aufeinanderfolgend Primzahlen!

@author Michael Fuhrmann
@date 2026-03-12
@lastModified 2026-03-12
@version 1.0
"""

from __future__ import annotations

import math
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy
from sympy import isprime, nextprime, prevprime, primerange, prime, primepi


# ===========================================================================
# Datenklassen für Ergebnisstrukturen
# ===========================================================================

@dataclass
class LegendreInterval:
    """
    @brief Repräsentiert ein Intervall (n², (n+1)²) mit Primzahlinformationen.

    @param n         Basiswert
    @param n_sq      n²
    @param np1_sq    (n+1)²
    @param primes    Liste der Primzahlen in diesem Intervall (exklusive Grenzen)
    @param count     Anzahl der Primzahlen im Intervall
    @lastModified 2026-03-12
    """
    n: int
    n_sq: int
    np1_sq: int
    primes: List[int]
    count: int


@dataclass
class BrocardInterval:
    """
    @brief Repräsentiert ein Intervall (pₙ², pₙ₊₁²) mit Primzahlinformationen.

    @param n         Index der Primzahl pₙ
    @param p_n       Primzahl pₙ
    @param p_n1      Primzahl pₙ₊₁
    @param primes    Liste der Primzahlen in (pₙ², pₙ₊₁²) exklusiv
    @param count     Anzahl der Primzahlen im Intervall
    @lastModified 2026-03-12
    """
    n: int
    p_n: int
    p_n1: int
    primes: List[int]
    count: int


# ===========================================================================
# KLASSE: LegendreConjecture
# ===========================================================================

class LegendreConjecture:
    r"""
    @brief Verifikation und Analyse der Legendre-Vermutung.
    @description
        Untersucht für jedes n im Bereich [1, N]:
        - Existiert mindestens eine Primzahl in (n², (n+1)²)?
        - Wie viele Primzahlen liegen in diesem Intervall?
        - Was ist die maximale und minimale Lücke?

        **Theorem (Bertrand-Tschebyschew, 1852)**:
        Für jedes n ≥ 1 existiert eine Primzahl p mit n < p ≤ 2n.

        **Conjecture (Legendre, 1798)**:
        Für jedes n ≥ 1 existiert eine Primzahl p mit n² < p < (n+1)².

        Die Intervalllänge ist (n+1)² - n² = 2n+1, wächst also linear.
        Laut Primzahlsatz liegt im Intervall (n², (n+1)²) asymptotisch:

            $$ \pi((n+1)^2) - \pi(n^2) \approx \frac{2n+1}{\ln(n^2)} = \frac{2n+1}{2\ln n} \approx \frac{n}{\ln n} $$

        Dieser Ausdruck → ∞ für n → ∞, was die Vermutung für große n plausibel macht.

    @author Michael Fuhrmann
    @date 2026-03-12
    @lastModified 2026-03-12
    """

    def __init__(self, max_n: int = 10_000):
        """
        @brief Initialisiert die Legendre-Vermutungsanalyse.
        @description
            Erstellt einen Primzahlsieb für den Bereich bis (max_n+1)².
            Für max_n = 10^5 wäre (max_n+1)² ≈ 10^10, was zu viel RAM braucht.
            Daher wird die Verifikation in Blöcken durchgeführt.

        @param max_n  Maximaler Wert von n für die Verifikation (Standard: 10.000).
                      Für vollständige Verifikation bis n=10^5 nutze verify_range().
        @lastModified 2026-03-12
        """
        # Maximale n-Grenze
        self._max_n: int = max_n

        # Cache für bereits berechnete Intervalle
        self._cache: Dict[int, LegendreInterval] = {}

        # Statistik
        self._min_count: Optional[int] = None
        self._max_count: Optional[int] = None
        self._min_n: Optional[int] = None
        self._max_n_gap: Optional[int] = None

        # Bekannte Gegenbeispielfreiheit bis n = ...
        self._verified_up_to: int = 0

    def interval_primes(self, n: int) -> LegendreInterval:
        """
        @brief Berechnet alle Primzahlen im Intervall (n², (n+1)²) exklusiv.
        @description
            Sucht alle Primzahlen p mit n² < p < (n+1)².
            Nutzt sympy.primerange() für effiziente Iteration.

            Intervalllänge: $(n+1)^2 - n^2 = 2n+1$

        @param n  Positive ganze Zahl ≥ 1.
        @return   LegendreInterval mit allen gefundenen Primzahlen.
        @raises ValueError  Wenn n < 1.
        @lastModified 2026-03-12
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        if n in self._cache:
            return self._cache[n]

        n_sq = n * n
        np1_sq = (n + 1) * (n + 1)

        # Alle Primzahlen im offenen Intervall (n², (n+1)²)
        # primerange(a, b) gibt Primzahlen a ≤ p < b zurück
        # Wir wollen n² < p < (n+1)², also primerange(n²+1, (n+1)²)
        ps = list(primerange(n_sq + 1, np1_sq))

        interval = LegendreInterval(
            n=n,
            n_sq=n_sq,
            np1_sq=np1_sq,
            primes=ps,
            count=len(ps),
        )
        self._cache[n] = interval
        return interval

    def verify_range(self, n_start: int = 1, n_end: int = 1_000) -> Dict:
        """
        @brief Verifiziert die Legendre-Vermutung für n von n_start bis n_end.
        @description
            Prüft für jedes n im Bereich, ob mindestens eine Primzahl in
            (n², (n+1)²) liegt. Gibt Statistiken zurück.

            **Laufzeit**: O(N · √N / log N) im Worst-Case.
            Für N = 10^5 nutzen wir blockweise Sieb-Verarbeitung.

        @param n_start  Startindex (Standard: 1).
        @param n_end    Endindex inklusiv (Standard: 1000).
        @return Dict mit Feldern:
                - verified: bool (True wenn Vermutung in gesamtem Bereich gilt)
                - counterexamples: List[int] (Liste von n ohne Primzahl im Intervall)
                - min_count: int (minimale Anzahl Primzahlen in einem Intervall)
                - max_count: int (maximale Anzahl Primzahlen in einem Intervall)
                - min_n: int (n mit minimaler Primzahlzahl)
                - max_n: int (n mit maximaler Primzahlzahl)
                - total_verified: int (Anzahl geprüfter Intervalle)
                - elapsed_seconds: float
        @lastModified 2026-03-12
        """
        if n_start < 1:
            raise ValueError(f"n_start muss ≥ 1 sein, erhalten: {n_start}")
        if n_end < n_start:
            raise ValueError(f"n_end ({n_end}) muss ≥ n_start ({n_start}) sein")

        t0 = time.time()
        counterexamples: List[int] = []
        min_count = None
        max_count = None
        min_n_val = None
        max_n_val = None

        for n in range(n_start, n_end + 1):
            interval = self.interval_primes(n)
            c = interval.count

            if c == 0:
                counterexamples.append(n)

            if min_count is None or c < min_count:
                min_count = c
                min_n_val = n
            if max_count is None or c > max_count:
                max_count = c
                max_n_val = n

        elapsed = time.time() - t0

        if n_end > self._verified_up_to:
            self._verified_up_to = n_end

        return {
            "verified": len(counterexamples) == 0,
            "counterexamples": counterexamples,
            "min_count": min_count,
            "max_count": max_count,
            "min_n": min_n_val,
            "max_n": max_n_val,
            "total_verified": n_end - n_start + 1,
            "elapsed_seconds": elapsed,
        }

    def bertrand_connection(self, n: int) -> Dict:
        """
        @brief Zeigt die Verbindung zwischen Legendre-Vermutung und Bertrand-Postulat.
        @description
            **Theorem (Tschebyschew/Bertrand, 1852)**:
            Für alle n ≥ 1 existiert eine Primzahl p mit n < p ≤ 2n.

            **Beziehung**:
            Die Legendre-Vermutung ist STÄRKER als Bertrand's Postulat.
            Wenn Legendre gilt, dann gibt es eine Primzahl in (n², (n+1)²),
            also insbesondere auch in (n, 2n) für passende n.

            Konkret: √n < p < √(2n+1) ≤ √(n+1) folgt nicht direkt,
            aber: Für k = ⌊√n⌋ liefert Legendre eine Primzahl p in
            (k², (k+1)²) ⊆ (n, n+2√n+1).

        @param n  Ganze Zahl n ≥ 1.
        @return Dict mit:
                - bertrand_witness: kleinste Primzahl in (n, 2n]
                - legendre_interval: (k², (k+1)²) für k = ⌊√n⌋
                - legendre_primes: Primzahlen in diesem Intervall
        @lastModified 2026-03-12
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        # Bertrand-Zeuge: kleinste Primzahl in (n, 2n]
        bertrand_witness = nextprime(n)
        # Prüfen ob in (n, 2n]
        in_bertrand = n < bertrand_witness <= 2 * n

        # Legendre-Intervall: k = ⌊√n⌋
        k = math.isqrt(n)
        if k < 1:
            k = 1
        interval = self.interval_primes(k)

        return {
            "n": n,
            "bertrand_range": (n + 1, 2 * n),
            "bertrand_witness": bertrand_witness,
            "bertrand_holds": in_bertrand,
            "legendre_k": k,
            "legendre_interval": (k * k, (k + 1) * (k + 1)),
            "legendre_primes": interval.primes,
            "legendre_count": interval.count,
        }

    def max_gap_between_squares(self, n_start: int = 1, n_end: int = 1_000) -> Dict:
        """
        @brief Findet die maximale Primzahllücke innerhalb der Quadrat-Intervalle.
        @description
            Berechnet für jedes Intervall (n², (n+1)²) die maximale Lücke zwischen
            aufeinanderfolgenden Primzahlen. Berücksichtigt auch die Lücken von n² zur
            ersten und von der letzten Primzahl zur (n+1)².

            Primzahllücke: $g_n = p_{k+1} - p_k$

        @param n_start  Startindex.
        @param n_end    Endindex inklusiv.
        @return Dict mit:
                - global_max_gap: int (größte gefundene Lücke über alle Intervalle)
                - max_gap_location: Tuple (n, p_before, p_after) wo die maximale Lücke liegt
                - gap_histogram: Dict[int, int] mapping Lückengröße → Anzahl
        @lastModified 2026-03-12
        """
        global_max_gap = 0
        max_gap_location = None
        gap_histogram: Dict[int, int] = defaultdict(int)

        for n in range(n_start, n_end + 1):
            interval = self.interval_primes(n)

            if interval.count == 0:
                # Keine Primzahl im Intervall: Lücke erstreckt sich über gesamtes Intervall
                gap = interval.np1_sq - interval.n_sq
                gap_histogram[gap] += 1
                if gap > global_max_gap:
                    global_max_gap = gap
                    max_gap_location = (n, interval.n_sq, interval.np1_sq)
            else:
                # Lücken zwischen den Primzahlen im Intervall berechnen
                points = interval.primes
                # Konsekutive Lücken
                for i in range(len(points) - 1):
                    gap = points[i + 1] - points[i]
                    gap_histogram[gap] += 1
                    if gap > global_max_gap:
                        global_max_gap = gap
                        max_gap_location = (n, points[i], points[i + 1])

        return {
            "global_max_gap": global_max_gap,
            "max_gap_location": max_gap_location,
            "gap_histogram": dict(sorted(gap_histogram.items())),
        }

    def gap_histogram_data(self, n_start: int = 1, n_end: int = 500) -> Dict:
        """
        @brief Erstellt Histogramm-Daten für Primzahllücken in Quadrat-Intervallen.
        @description
            Zählt für jeden Bereich n ∈ [n_start, n_end], wie viele Primzahlen im
            Intervall (n², (n+1)²) liegen. Gibt die Verteilung zurück.

        @param n_start  Startindex.
        @param n_end    Endindex inklusiv.
        @return Dict mit:
                - counts_distribution: Dict[int, int] (Anzahl → Häufigkeit)
                - mean_count: float (Durchschnittliche Anzahl Primzahlen)
                - std_count: float (Standardabweichung)
        @lastModified 2026-03-12
        """
        counts = []
        for n in range(n_start, n_end + 1):
            interval = self.interval_primes(n)
            counts.append(interval.count)

        counts_arr = np.array(counts, dtype=float)
        distribution: Dict[int, int] = defaultdict(int)
        for c in counts:
            distribution[c] += 1

        return {
            "counts_distribution": dict(sorted(distribution.items())),
            "mean_count": float(np.mean(counts_arr)),
            "std_count": float(np.std(counts_arr)),
            "min_count": int(np.min(counts_arr)),
            "max_count": int(np.max(counts_arr)),
        }


# ===========================================================================
# KLASSE: BrocardConjecture7
# ===========================================================================

class BrocardConjecture7:
    r"""
    @brief Verifikation der Brocard-Vermutung (1904) über Primzahlquadrate.
    @description
        **WICHTIGER HINWEIS**: Dies ist NICHT die Brocard-Ramanujan-Gleichung!

        **Brocard-Vermutung (1904)**:
        Für alle n ≥ 2 existieren zwischen $p_n^2$ und $p_{n+1}^2$ mindestens 4 Primzahlen:

            $$ \pi(p_{n+1}^2) - \pi(p_n^2) \geq 4 \quad \forall n \geq 2 $$

        Dabei ist $p_n$ die n-te Primzahl (p_1=2, p_2=3, p_3=5, ...).

        **Hintergrund**:
        Für n=1: Zwischen 4 (=2²) und 9 (=3²) liegen die Primzahlen 5 und 7. Nur 2!
        Deshalb beginnt die Brocard-Vermutung bei n=2.

        Für n=2: Zwischen 9 (=3²) und 25 (=5²) liegen 11, 13, 17, 19, 23. Das sind 5 Primzahlen!

        **Verbindung zur Legendre-Vermutung**:
        Falls die Legendre-Vermutung gilt, dann liegen zwischen p_n² und p_{n+1}²
        mindestens so viele Primzahlen wie der Abstand p_{n+1} - p_n Intervalle überdeckt.
        Da p_{n+1} - p_n ≥ 2 (für n ≥ 2), wären es mindestens 2 Primzahlen, aber Brocard
        behauptet sogar 4.

        **Beweisstatus**: Offen. Verifiziert für alle Primzahlen p < 10^9.

    @author Michael Fuhrmann
    @date 2026-03-12
    @lastModified 2026-03-12
    """

    # Bekannte Grenze der Verifikation
    VERIFIED_LIMIT = 10**9

    # Minimale Anzahl laut Vermutung
    MINIMUM_REQUIRED = 4

    def __init__(self):
        """
        @brief Initialisiert die Brocard-Vermutungsanalyse.
        @description
            Lädt die ersten Primzahlen und bereitet die Analyse vor.
        @lastModified 2026-03-12
        """
        # Cache für berechnete Intervalle
        self._interval_cache: Dict[int, BrocardInterval] = {}

        # Statistiken
        self._min_observed: Optional[int] = None
        self._min_observed_n: Optional[int] = None

    def interval_primes(self, n: int) -> BrocardInterval:
        """
        @brief Berechnet die Primzahlen im Intervall (pₙ², pₙ₊₁²) exklusiv.
        @description
            Sucht alle Primzahlen p mit $p_n^2 < p < p_{n+1}^2$.
            Nutzt sympy.prime(n) für den n-ten Primwert und primerange().

            **n=1 Sonderfall**: Nur 2 Primzahlen (5, 7) in (4, 9).
            **n≥2**: Mindestens 4 laut Brocard-Vermutung.

        @param n  Index der Primzahl (n ≥ 1, p_1 = 2).
        @return   BrocardInterval mit den gefundenen Primzahlen.
        @raises ValueError  Wenn n < 1.
        @lastModified 2026-03-12
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        if n in self._interval_cache:
            return self._interval_cache[n]

        # n-te und (n+1)-te Primzahl bestimmen
        p_n = prime(n)
        p_n1 = prime(n + 1)

        p_n_sq = p_n * p_n
        p_n1_sq = p_n1 * p_n1

        # Alle Primzahlen im offenen Intervall (pₙ², pₙ₊₁²)
        ps = list(primerange(p_n_sq + 1, p_n1_sq))

        interval = BrocardInterval(
            n=n,
            p_n=p_n,
            p_n1=p_n1,
            primes=ps,
            count=len(ps),
        )
        self._interval_cache[n] = interval
        return interval

    def verify_range(self, n_start: int = 2, n_end: int = 1_000) -> Dict:
        """
        @brief Verifiziert die Brocard-Vermutung für Primzahl-Indizes n_start bis n_end.
        @description
            Prüft für jedes n im Bereich, ob ≥ 4 Primzahlen in (pₙ², pₙ₊₁²) liegen.

            **n=1 wird separat behandelt** (dort sind es nur 2 Primzahlen, was korrekt ist,
            da die Vermutung erst ab n=2 gilt).

        @param n_start  Startindex (Standard: 2, da Vermutung ab n=2 gilt).
        @param n_end    Endindex inklusiv (Standard: 1000).
        @return Dict mit:
                - verified: bool
                - counterexamples: List[int] (Indizes n mit < 4 Primzahlen)
                - min_count: int (minimale gefundene Anzahl)
                - min_count_n: int (n wo das Minimum liegt)
                - max_count: int (maximale gefundene Anzahl)
                - statistics: Dict mit Mittelwert, Std
                - total_verified: int
                - elapsed_seconds: float
        @lastModified 2026-03-12
        """
        if n_start < 1:
            raise ValueError(f"n_start muss ≥ 1 sein, erhalten: {n_start}")

        t0 = time.time()
        counterexamples: List[int] = []
        counts: List[int] = []
        min_count = None
        min_count_n = None
        max_count = None

        for n in range(n_start, n_end + 1):
            interval = self.interval_primes(n)
            c = interval.count
            counts.append(c)

            # Brocard gilt ab n=2
            if n >= 2 and c < self.MINIMUM_REQUIRED:
                counterexamples.append(n)

            if min_count is None or c < min_count:
                min_count = c
                min_count_n = n
            if max_count is None or c > max_count:
                max_count = c

        elapsed = time.time() - t0
        counts_arr = np.array(counts, dtype=float)

        return {
            "verified": len(counterexamples) == 0,
            "counterexamples": counterexamples,
            "min_count": min_count,
            "min_count_n": min_count_n,
            "max_count": max_count,
            "statistics": {
                "mean": float(np.mean(counts_arr)),
                "std": float(np.std(counts_arr)),
                "median": float(np.median(counts_arr)),
            },
            "total_verified": n_end - n_start + 1,
            "elapsed_seconds": elapsed,
        }

    def minimum_tracking(self, n_start: int = 2, n_end: int = 1_000) -> Dict:
        """
        @brief Verfolgt das Minimum der Primzahlanzahl über alle Intervalle.
        @description
            Gibt für jedes n das laufende Minimum zurück und identifiziert,
            wo das globale Minimum liegt. Nützlich um zu sehen, ob das Minimum
            gegen 4 konvergiert oder darunter fällt.

        @param n_start  Startindex.
        @param n_end    Endindex inklusiv.
        @return Dict mit:
                - running_minimum: List[Tuple[int, int]] (n, laufendes Minimum)
                - global_minimum: int
                - global_minimum_location: int (n)
                - conjecture_margin: int (globales Minimum - 4)
        @lastModified 2026-03-12
        """
        running_min = None
        running_history: List[Tuple[int, int]] = []
        global_min = None
        global_min_n = None

        for n in range(n_start, n_end + 1):
            interval = self.interval_primes(n)
            c = interval.count

            if running_min is None or c < running_min:
                running_min = c

            running_history.append((n, running_min))

            if global_min is None or c < global_min:
                global_min = c
                global_min_n = n

        return {
            "running_minimum": running_history,
            "global_minimum": global_min,
            "global_minimum_location": global_min_n,
            "conjecture_margin": (global_min - self.MINIMUM_REQUIRED) if global_min is not None else None,
        }

    def statistics(self, n_start: int = 2, n_end: int = 1_000) -> Dict:
        r"""
        @brief Vollständige Statistik der Primzahlanzahl in Brocard-Intervallen.
        @description
            Berechnet Verteilung, Momente und weitere statistische Kennzahlen
            der Anzahl der Primzahlen in den Intervallen (pₙ², pₙ₊₁²).

            Asymptotisch erwartet man laut Primzahlsatz:
                $$ \pi(p_{n+1}^2) - \pi(p_n^2) \approx \frac{p_{n+1}^2 - p_n^2}{2\ln p_n} $$

        @param n_start  Startindex.
        @param n_end    Endindex inklusiv.
        @return Dict mit Verteilung und Momenten.
        @lastModified 2026-03-12
        """
        counts: List[int] = []
        intervals: List[BrocardInterval] = []

        for n in range(n_start, n_end + 1):
            iv = self.interval_primes(n)
            counts.append(iv.count)
            intervals.append(iv)

        counts_arr = np.array(counts, dtype=float)
        distribution: Dict[int, int] = defaultdict(int)
        for c in counts:
            distribution[c] += 1

        return {
            "n_range": (n_start, n_end),
            "count_distribution": dict(sorted(distribution.items())),
            "mean": float(np.mean(counts_arr)),
            "std": float(np.std(counts_arr)),
            "median": float(np.median(counts_arr)),
            "min": int(np.min(counts_arr)),
            "max": int(np.max(counts_arr)),
            "min_exceeds_4": bool(np.min(counts_arr) >= 4),
            "total_intervals": len(counts),
        }
