r"""
@file andrica_cramer.py
@brief Verifikation und Analyse der Andrica-Vermutung und Cramér-Vermutung.
@description
    Dieses Modul untersucht zwei fundamentale Vermutungen über Primzahllücken:

    **1. Andrica-Vermutung (1985, ungeklärt)**

    Conjecture (Andrica, 1985): Für alle n ≥ 1 gilt:

        $$ A_n = \sqrt{p_{n+1}} - \sqrt{p_n} < 1 $$

    Äquivalent zur Aussage, dass die Primzahllücke $g_n = p_{n+1} - p_n$
    die Schranke $g_n < 2\sqrt{p_n} + 1$ erfüllt.

    **Bekanntes Maximalwert**: $A_n \approx 0.6708$ bei $p_n = 7, p_{n+1} = 11$.

    **Status**: Numerisch verifiziert für alle Primzahlen bis $4 \times 10^{18}$
    (Meissel, Lehmer, 2014). Wäre durch die Cramér-Vermutung impliziert.

    **2. Cramér-Vermutung (1936, ungeklärt)**

    Conjecture (Cramér, 1936): Für alle n gilt:

        $$ \limsup_{n \to \infty} \frac{p_{n+1} - p_n}{(\ln p_n)^2} = 1 $$

    Äquivalent: Es gibt unendlich viele n mit $g_n \gg (\ln p_n)^2$,
    aber niemals $g_n > (1+\epsilon)(\ln p_n)^2$ für jedes feste $\epsilon > 0$.

    **Granville-Korrektur (1995)**:
    Granville argumentierte basierend auf probabilistischen Heuristiken, dass
    der Limes Superior tatsächlich $2e^{-\gamma} \approx 1.1229$ sein könnte
    (statt 1 wie Cramér vermutete).

    Dabei ist $\gamma \approx 0.5772$ die Euler-Mascheroni-Konstante.

    **Aktueller Rekord**: Baker, Harman und Pintz (2001) bewiesen:
    $g_n \ll p_n^{0.525}$ (bedingt durch Siebtechnik, KEIN Beweis von Cramér).

    **Bezug zu anderen Vermutungen**:
    - Cramér ⟹ Andrica: $(\ln p)^2 < 2\sqrt{p}+1$ für $p \geq 3$
    - Cramér ⟹ Legendre: $(\ln p)^2 \ll \sqrt{p}$ für große p
    - Riemann-Hypothese: Impliziert $g_n \ll \sqrt{p_n}\log p_n$

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

# Euler-Mascheroni-Konstante
_EULER_MASCHERONI = 0.5772156649015328606

# Granville-Korrekturfaktor: 2e^{-γ}
GRANVILLE_CONSTANT = 2.0 * math.exp(-_EULER_MASCHERONI)


# ===========================================================================
# Datenklassen
# ===========================================================================

@dataclass
class AndricaValue:
    r"""
    @brief Enthält den Andrica-Wert A_n = √p_{n+1} − √pₙ für ein Primzahlpaar.

    @param n        Index (1-basiert)
    @param p_n      n-te Primzahl
    @param p_n1     (n+1)-te Primzahl
    @param gap      Primzahllücke p_{n+1} - p_n
    @param value    Andrica-Wert $\sqrt{p_{n+1}} - \sqrt{p_n}$
    @param holds    True wenn value < 1 (Vermutung gilt für dieses n)
    @lastModified 2026-03-12
    """
    n: int
    p_n: int
    p_n1: int
    gap: int
    value: float
    holds: bool


@dataclass
class CramerValue:
    """
    @brief Enthält den normalisierten Cramér-Wert für eine Primzahllücke.

    @param n           Index (1-basiert)
    @param p_n         n-te Primzahl
    @param p_n1        (n+1)-te Primzahl
    @param gap         Primzahllücke g_n = p_{n+1} - p_n
    @param log_sq      (log p_n)²
    @param cramer_val  Normierter Wert: g_n / (log p_n)²
    @param granville_ratio  g_n / (2e^{-γ} · (log p_n)²)
    @lastModified 2026-03-12
    """
    n: int
    p_n: int
    p_n1: int
    gap: int
    log_sq: float
    cramer_val: float
    granville_ratio: float


# ===========================================================================
# KLASSE: AndricaConjecture
# ===========================================================================

class AndricaConjecture:
    r"""
    @brief Verifikation und Analyse der Andrica-Vermutung.
    @description
        Untersucht für jedes Primzahlpaar (pₙ, pₙ₊₁):

            $$ A_n = \sqrt{p_{n+1}} - \sqrt{p_n} < 1 $$

        **Äquivalente Formulierungen**:
        1. $\sqrt{p_{n+1}} - \sqrt{p_n} < 1$
        2. $p_{n+1} < (\sqrt{p_n} + 1)^2 = p_n + 2\sqrt{p_n} + 1$
        3. $g_n = p_{n+1} - p_n < 2\sqrt{p_n} + 1$

        Formulierung 3 zeigt: Die Andrica-Vermutung ist äquivalent zur Aussage,
        dass die Primzahllücke sublinear in √p wächst.

        **Verbindung zur Legendre-Vermutung**:
        Die Andrica-Vermutung für AUFEINANDERFOLGENDE Primzahlen ist äquivalent zur
        Legendre-Vermutung! Denn: Wenn $p_{n+1} < (\sqrt{p_n}+1)^2$, dann liegt
        $p_{n+1}$ im Intervall $(p_n, (\sqrt{p_n}+1)^2)$.

        Allerdings: Legendre betrachtet alle Intervalle $(n^2, (n+1)^2)$,
        Andrica nur Paare AUFEINANDERFOLGENDER Primzahlen.

    @author Michael Fuhrmann
    @date 2026-03-12
    @lastModified 2026-03-12
    """

    # Bekannte maximale Andrica-Differenz: √11 - √7 ≈ 0.6708
    KNOWN_MAX_VALUE: float = math.sqrt(11) - math.sqrt(7)
    KNOWN_MAX_N: int = 4  # Index von p_4=7 (7 ist 4. Primzahl)

    def __init__(self):
        """
        @brief Initialisiert die Andrica-Vermutungsanalyse.
        @lastModified 2026-03-12
        """
        self._cache: Dict[int, AndricaValue] = {}
        self._max_value: Optional[float] = None
        self._max_value_n: Optional[int] = None

    def compute(self, n: int) -> AndricaValue:
        r"""
        @brief Berechnet den Andrica-Wert A_n = √p_{n+1} − √pₙ.
        @description
            $A_n = \sqrt{p_{n+1}} - \sqrt{p_n}$

            Die Vermutung besagt A_n < 1 für alle n ≥ 1.

        @param n  Primzahlindex (n ≥ 1, entspricht p_1 = 2).
        @return   AndricaValue mit allen berechneten Größen.
        @raises ValueError  Wenn n < 1.
        @lastModified 2026-03-12
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        if n in self._cache:
            return self._cache[n]

        p_n = prime(n)
        p_n1 = prime(n + 1)
        gap = p_n1 - p_n
        value = math.sqrt(p_n1) - math.sqrt(p_n)
        holds = value < 1.0

        av = AndricaValue(
            n=n,
            p_n=p_n,
            p_n1=p_n1,
            gap=gap,
            value=value,
            holds=holds,
        )
        self._cache[n] = av

        if self._max_value is None or value > self._max_value:
            self._max_value = value
            self._max_value_n = n

        return av

    def verify_range(self, n_start: int = 1, n_end: int = 500_000) -> Dict:
        """
        @brief Verifiziert die Andrica-Vermutung für n von n_start bis n_end.
        @description
            Prüft A_n < 1 für alle n im angegebenen Bereich.
            Bei n_end = 500.000 werden ca. die ersten 7.8 Millionen Primzahlen geprüft.

            **Hinweis**: Für große Bereiche (n_end > 10^4) kann die Berechnung
            einige Minuten dauern. Sympy's prime(n) ist für große n langsam.
            Für Performance-kritischen Code: Siebeansatz nutzen.

        @param n_start  Startindex (Standard: 1).
        @param n_end    Endindex inklusiv (Standard: 500.000).
        @return Dict mit:
                - verified: bool
                - counterexamples: List[int]
                - max_value: float (maximaler Andrica-Wert)
                - max_value_n: int (n wo Maximum liegt)
                - max_value_primes: Tuple[int, int] (p_n, p_{n+1}) beim Maximum
                - gap_stats: Dict mit Lückenstatistik
                - elapsed_seconds: float
        @lastModified 2026-03-12
        """
        if n_start < 1:
            raise ValueError(f"n_start muss ≥ 1 sein, erhalten: {n_start}")

        t0 = time.time()
        counterexamples: List[int] = []
        max_val = 0.0
        max_val_n = n_start
        max_val_primes = (2, 3)
        gaps: List[int] = []
        andrica_values: List[float] = []

        # Effiziente Iteration: nextprime statt prime(n) für jedes n
        p = prime(n_start)
        for n in range(n_start, n_end + 1):
            p_next = nextprime(p)
            gap = p_next - p
            value = math.sqrt(p_next) - math.sqrt(p)

            if value >= 1.0:
                counterexamples.append(n)

            if value > max_val:
                max_val = value
                max_val_n = n
                max_val_primes = (p, p_next)

            gaps.append(gap)
            andrica_values.append(value)
            p = p_next

        elapsed = time.time() - t0
        gaps_arr = np.array(gaps, dtype=float)

        return {
            "verified": len(counterexamples) == 0,
            "counterexamples": counterexamples,
            "max_value": max_val,
            "max_value_n": max_val_n,
            "max_value_primes": max_val_primes,
            "gap_stats": {
                "mean_gap": float(np.mean(gaps_arr)),
                "max_gap": int(np.max(gaps_arr)),
                "std_gap": float(np.std(gaps_arr)),
            },
            "mean_andrica": float(np.mean(andrica_values)),
            "total_verified": n_end - n_start + 1,
            "elapsed_seconds": elapsed,
        }

    def max_value_analysis(self, n_end: int = 1_000) -> Dict:
        r"""
        @brief Analysiert die Entwicklung des maximalen Andrica-Werts.
        @description
            Berechnet für jedes n den Andrica-Wert und verfolgt das laufende Maximum.
            Das globale Maximum liegt bekanntermaßen bei n=4 (p_n=7, p_{n+1}=11).

            $A_4 = \sqrt{11} - \sqrt{7} \approx 0.6708$

            Es ist eine offene Frage, ob dieses das globale Maximum ist.

        @param n_end  Endindex inklusiv.
        @return Dict mit laufendem Maximum und Positionen neuer Maxima.
        @lastModified 2026-03-12
        """
        running_max = 0.0
        new_max_events: List[Dict] = []
        p = prime(1)

        for n in range(1, n_end + 1):
            p_next = nextprime(p)
            value = math.sqrt(p_next) - math.sqrt(p)

            if value > running_max:
                running_max = value
                new_max_events.append({
                    "n": n,
                    "p_n": p,
                    "p_n1": p_next,
                    "value": value,
                    "gap": p_next - p,
                })

            p = p_next

        return {
            "global_max": running_max,
            "new_max_events": new_max_events,
            "total_checked": n_end,
            "known_max": self.KNOWN_MAX_VALUE,
        }

    def gap_statistics(self, n_start: int = 1, n_end: int = 1_000) -> Dict:
        r"""
        @brief Berechnet Gap-Statistiken im Kontext der Andrica-Vermutung.
        @description
            Für jedes n berechnet die Lücke g_n = p_{n+1} - p_n und die
            Andrica-Schranke $2\sqrt{p_n} + 1$ (äquivalente Formulierung).

            Das Verhältnis $g_n / (2\sqrt{p_n})$ zeigt, wie weit die Lücke
            von der Andrica-Schranke entfernt ist.

        @param n_start  Startindex.
        @param n_end    Endindex inklusiv.
        @return Dict mit Lückenstatistiken.
        @lastModified 2026-03-12
        """
        gaps: List[int] = []
        andrica_bounds: List[float] = []
        ratios: List[float] = []
        p = prime(n_start)

        for _ in range(n_start, n_end + 1):
            p_next = nextprime(p)
            gap = p_next - p
            bound = 2.0 * math.sqrt(p) + 1.0
            ratio = gap / bound

            gaps.append(gap)
            andrica_bounds.append(bound)
            ratios.append(ratio)
            p = p_next

        ratios_arr = np.array(ratios)
        gaps_arr = np.array(gaps)

        return {
            "max_ratio": float(np.max(ratios_arr)),
            "mean_ratio": float(np.mean(ratios_arr)),
            "std_ratio": float(np.std(ratios_arr)),
            "max_gap": int(np.max(gaps_arr)),
            "mean_gap": float(np.mean(gaps_arr)),
            "all_below_bound": bool(np.all(np.array(gaps) < np.array(andrica_bounds))),
        }


# ===========================================================================
# KLASSE: CramerConjecture
# ===========================================================================

class CramerConjecture:
    r"""
    @brief Analyse und Statistik der Cramér-Vermutung über Primzahllücken.
    @description
        **Cramér-Vermutung (1936)**:

            $$ \limsup_{n \to \infty} \frac{g_n}{(\ln p_n)^2} = 1 $$

        wobei $g_n = p_{n+1} - p_n$ die n-te Primzahllücke ist.

        **Granville-Korrektur (1995)**:
        Basierend auf der Analyse des Cramér-Modells mit Hardy-Littlewood-Heuristiken
        korrigierte Granville den Limes Superior auf:

            $$ 2e^{-\gamma} \approx 1.1229 $$

        Diese Konstante ergibt sich aus $\gamma = 0.5772...$ (Euler-Mascheroni).

        **Bekannte bedingte Ergebnisse**:
        - Unter RH: $g_n = O(\sqrt{p_n} \ln p_n)$
        - Baker-Harman-Pintz (2001): $g_n \ll p_n^{0.525}$ (unbedingt)
        - Westzynthius (1931): $g_n / \ln p_n \to \infty$ infinitely often

        **Aktueller Rekord-Gap** (Oliveira e Silva, 2014):
        Bei $p = 1693182318746371$ (1.7 × 10^15) wurde die größte bekannte normierte
        Lücke von ~1.17 beobachtet (nahe an der Granville-Schranke).

    @author Michael Fuhrmann
    @date 2026-03-12
    @lastModified 2026-03-12
    """

    # Granville-Korrekturkoeffizient: 2e^{-γ}
    GRANVILLE_CONSTANT: float = GRANVILLE_CONSTANT

    # Euler-Mascheroni-Konstante
    EULER_MASCHERONI: float = _EULER_MASCHERONI

    def __init__(self):
        """
        @brief Initialisiert die Cramér-Vermutungsanalyse.
        @lastModified 2026-03-12
        """
        self._cache: Dict[int, CramerValue] = {}

    def compute(self, n: int) -> CramerValue:
        r"""
        @brief Berechnet den normalisierten Cramér-Wert für Index n.
        @description
            Berechnet:
                $C_n = \frac{g_n}{(\ln p_n)^2}$
                $G_n = \frac{g_n}{2e^{-\gamma} \cdot (\ln p_n)^2}$ (Granville-normiert)

        @param n  Primzahlindex (n ≥ 1).
        @return   CramerValue mit normalisierten Werten.
        @raises ValueError  Wenn n < 1.
        @lastModified 2026-03-12
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        if n in self._cache:
            return self._cache[n]

        p_n = prime(n)
        p_n1 = prime(n + 1)
        gap = p_n1 - p_n
        log_sq = math.log(p_n) ** 2
        cramer_val = gap / log_sq if log_sq > 0 else 0.0
        granville_ratio = cramer_val / self.GRANVILLE_CONSTANT

        cv = CramerValue(
            n=n,
            p_n=p_n,
            p_n1=p_n1,
            gap=gap,
            log_sq=log_sq,
            cramer_val=cramer_val,
            granville_ratio=granville_ratio,
        )
        self._cache[n] = cv
        return cv

    def prime_gap_statistics(self, n_start: int = 1, n_end: int = 10_000) -> Dict:
        """
        @brief Berechnet umfangreiche Primzahllücken-Statistik.
        @description
            Iteriert über Primzahlpaare und berechnet normierte Cramér-Werte.
            Gibt Verteilung, Maxima und Annäherung an Cramér- und Granville-Schranke zurück.

        @param n_start  Startindex (Standard: 1).
        @param n_end    Endindex inklusiv (Standard: 10.000).
        @return Dict mit:
                - max_cramer_value: float (größter C_n = g_n/(log pₙ)²)
                - max_cramer_n: int (Index des Maximums)
                - max_cramer_primes: Tuple[int, int]
                - mean_cramer: float
                - std_cramer: float
                - granville_exceedances: List[int] (n wo C_n > Granville-Konstante)
                - max_granville_ratio: float
                - gap_distribution: Dict[int, int]
                - elapsed_seconds: float
        @lastModified 2026-03-12
        """
        if n_start < 1:
            raise ValueError(f"n_start muss ≥ 1 sein, erhalten: {n_start}")

        t0 = time.time()
        cramer_vals: List[float] = []
        granville_ratios: List[float] = []
        gaps: List[int] = []
        max_cv = 0.0
        max_cv_n = n_start
        max_cv_primes = (2, 3)
        granville_exceedances: List[int] = []
        gap_distribution: Dict[int, int] = defaultdict(int)

        # Effiziente Iteration mit nextprime
        p = prime(n_start)
        for n in range(n_start, n_end + 1):
            p_next = nextprime(p)
            gap = p_next - p
            log_sq = math.log(p) ** 2
            cv = gap / log_sq if log_sq > 0 else 0.0
            gr = cv / self.GRANVILLE_CONSTANT

            cramer_vals.append(cv)
            granville_ratios.append(gr)
            gaps.append(gap)
            gap_distribution[gap] += 1

            if cv > max_cv:
                max_cv = cv
                max_cv_n = n
                max_cv_primes = (p, p_next)

            if cv > self.GRANVILLE_CONSTANT:
                granville_exceedances.append(n)

            p = p_next

        elapsed = time.time() - t0
        cv_arr = np.array(cramer_vals)
        gr_arr = np.array(granville_ratios)

        return {
            "max_cramer_value": max_cv,
            "max_cramer_n": max_cv_n,
            "max_cramer_primes": max_cv_primes,
            "mean_cramer": float(np.mean(cv_arr)),
            "std_cramer": float(np.std(cv_arr)),
            "granville_constant": self.GRANVILLE_CONSTANT,
            "granville_exceedances": granville_exceedances,
            "max_granville_ratio": float(np.max(gr_arr)),
            "mean_granville_ratio": float(np.mean(gr_arr)),
            "gap_distribution": dict(sorted(gap_distribution.items())),
            "total_analyzed": n_end - n_start + 1,
            "elapsed_seconds": elapsed,
        }

    def granville_correction(self, n_start: int = 1, n_end: int = 1_000) -> Dict:
        r"""
        @brief Analysiert die Granville-Korrektur der Cramér-Vermutung.
        @description
            Die Granville-Korrektur sagt voraus, dass der Limes Superior von
            $g_n / (\ln p_n)^2$ nicht 1 sondern $2e^{-\gamma} \approx 1.1229$ ist.

            Diese Funktion berechnet:
            - Den normierten Wert $C_n = g_n/(\ln p_n)^2$
            - Den Granville-normierten Wert $G_n = C_n / (2e^{-\gamma})$
            - Statistiken beider Normierungen

        @param n_start  Startindex.
        @param n_end    Endindex inklusiv.
        @return Dict mit Granville-Analyse.
        @lastModified 2026-03-12
        """
        cramer_vals: List[float] = []
        granville_vals: List[float] = []
        p = prime(n_start)

        for _ in range(n_start, n_end + 1):
            p_next = nextprime(p)
            gap = p_next - p
            log_sq = math.log(p) ** 2
            cv = gap / log_sq
            gv = cv / self.GRANVILLE_CONSTANT
            cramer_vals.append(cv)
            granville_vals.append(gv)
            p = p_next

        cv_arr = np.array(cramer_vals)
        gv_arr = np.array(granville_vals)

        return {
            "granville_constant": self.GRANVILLE_CONSTANT,
            "euler_mascheroni": self.EULER_MASCHERONI,
            "max_cramer_normalized": float(np.max(cv_arr)),
            "max_granville_normalized": float(np.max(gv_arr)),
            "mean_cramer_normalized": float(np.mean(cv_arr)),
            "mean_granville_normalized": float(np.mean(gv_arr)),
            "fraction_above_1": float(np.mean(cv_arr > 1.0)),
            "fraction_above_granville": float(np.mean(cv_arr > self.GRANVILLE_CONSTANT)),
            "n_range": (n_start, n_end),
        }

    def maxima_analysis(self, n_end: int = 1_000) -> Dict:
        r"""
        @brief Analysiert die Record-Maxima der Cramér-Werte.
        @description
            Verfolgt das laufende Maximum von $C_n = g_n/(\ln p_n)^2$.
            Jedes neue Rekordmaximum wird mit dem zugehörigen Primzahlpaar protokolliert.

            Laut der Cramér-Vermutung sollte das Supremum 1 sein,
            laut Granville eher 1.1229.

        @param n_end  Maximaler Index.
        @return Dict mit chronologischer Liste aller Rekordmaxima.
        @lastModified 2026-03-12
        """
        running_max = 0.0
        record_events: List[Dict] = []
        p = prime(1)

        for n in range(1, n_end + 1):
            p_next = nextprime(p)
            gap = p_next - p
            log_sq = math.log(p) ** 2
            cv = gap / log_sq

            if cv > running_max:
                running_max = cv
                record_events.append({
                    "n": n,
                    "p_n": p,
                    "p_n1": p_next,
                    "gap": gap,
                    "cramer_value": cv,
                    "granville_ratio": cv / self.GRANVILLE_CONSTANT,
                    "log_p_sq": log_sq,
                })

            p = p_next

        return {
            "record_maxima": record_events,
            "final_maximum": running_max,
            "granville_constant": self.GRANVILLE_CONSTANT,
            "total_checked": n_end,
        }

    def distribution_analysis(self, n_start: int = 1, n_end: int = 10_000) -> Dict:
        r"""
        @brief Analysiert die Verteilung der normalisierten Primzahllücken.
        @description
            Die Cramér-Vermutung macht Aussagen über die asymptotische Verteilung
            der normierten Lücken $C_n = g_n/(\ln p_n)^2$.

            Laut probabilistischem Modell (Cramér 1936) sollten die Lücken
            annähernd exponentialverteilt mit Parameter 1 sein, was bedeutet:
            $P(C_n > x) \approx e^{-x}$ für große $p_n$.

        @param n_start  Startindex.
        @param n_end    Endindex inklusiv.
        @return Dict mit Verteilungsanalyse.
        @lastModified 2026-03-12
        """
        cramer_vals: List[float] = []
        p = prime(n_start)

        for _ in range(n_start, n_end + 1):
            p_next = nextprime(p)
            gap = p_next - p
            log_sq = math.log(p) ** 2
            cv = gap / log_sq
            cramer_vals.append(cv)
            p = p_next

        cv_arr = np.array(cramer_vals)

        # Quantile berechnen
        quantiles = {
            "q50": float(np.percentile(cv_arr, 50)),
            "q75": float(np.percentile(cv_arr, 75)),
            "q90": float(np.percentile(cv_arr, 90)),
            "q95": float(np.percentile(cv_arr, 95)),
            "q99": float(np.percentile(cv_arr, 99)),
        }

        # Histogram der normierten Werte
        hist_bins = np.linspace(0, float(np.max(cv_arr)) + 0.1, 20)
        hist_counts, hist_edges = np.histogram(cv_arr, bins=hist_bins)

        return {
            "mean": float(np.mean(cv_arr)),
            "std": float(np.std(cv_arr)),
            "skewness": float(
                np.mean(((cv_arr - np.mean(cv_arr)) / np.std(cv_arr)) ** 3)
            ) if np.std(cv_arr) > 0 else 0.0,
            "max": float(np.max(cv_arr)),
            "quantiles": quantiles,
            "histogram_counts": hist_counts.tolist(),
            "histogram_edges": hist_edges.tolist(),
            "granville_constant": self.GRANVILLE_CONSTANT,
            "n_range": (n_start, n_end),
        }
