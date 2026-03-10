"""
@file additive_number_theory.py
@brief Additive Zahlentheorie: Partitionen, Waring-Problem, Goldbach, AP, Summenmengen.

@description
    Dieses Modul implementiert die wichtigsten Gebiete der additiven Zahlentheorie:

    1. PartitionFunction       – Partitionen ganzer Zahlen, Euler-Pentagonaltheorem,
                                  Ramanujan-Kongruenzen
    2. HardyRamanujanCircleMethod – Asymptotik und Rademacher-Formel für p(n)
    3. WaringProblem           – Darstellung als Summe k-ter Potenzen (Lagrange, Legendre)
    4. GoldbachAnalysis        – Goldbach-Zerlegungen, Komet-Daten, schwache Vermutung
    5. ArithmeticProgressions  – Dirichlet, Van-der-Waerden, Green-Tao
    6. SumSetTheory            – Summenmengen, Cauchy-Davenport, additive Energie
    7. Standalone-Funktionen   – Schur-Zahlen, Happy-Ending, Bertrand-Postulat

    Alle Algorithmen sind ausführlich kommentiert und dienen gleichzeitig als
    Lernmaterial für additive Kombinatorik und analytische Zahlentheorie.

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import math
import itertools
import functools
from typing import Generator, Optional
from fractions import Fraction


# ===========================================================================
# HILFSFUNKTIONEN (intern)
# ===========================================================================

def _sieve_primes(limit: int) -> list[int]:
    """
    Sieb des Eratosthenes bis zur angegebenen Obergrenze.

    @param limit: Obergrenze (inklusiv)
    @return: Liste aller Primzahlen ≤ limit
    @lastModified 2026-03-10
    """
    if limit < 2:
        return []
    # Boolesches Array: is_p[i] == True ⟺ i ist prim
    is_p = bytearray([1]) * (limit + 1)
    is_p[0] = is_p[1] = 0
    p = 2
    while p * p <= limit:
        if is_p[p]:
            # Alle Vielfachen von p als zusammengesetzt markieren
            is_p[p * p::p] = bytearray(len(is_p[p * p::p]))
        p += 1
    return [i for i, v in enumerate(is_p) if v]


def _is_prime_simple(n: int) -> bool:
    """
    Einfacher Primalitätstest via Probedivision.

    @param n: Zu testende Zahl
    @return: True wenn n prim
    @lastModified 2026-03-10
    """
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    # Probedivision mit 6k±1-Schritten
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def _gcd(a: int, b: int) -> int:
    """
    Größter gemeinsamer Teiler via Euklidischem Algorithmus.

    @param a: Erste Zahl
    @param b: Zweite Zahl
    @return: ggT(a, b)
    @lastModified 2026-03-10
    """
    while b:
        a, b = b, a % b
    return abs(a)


# ===========================================================================
# KLASSE 1: PARTITIONSFUNKTION
# ===========================================================================

class PartitionFunction:
    """
    Klasse für die Partitionsfunktion p(n) und verwandte Konzepte.

    Eine Partition von n ist eine Zerlegung von n in positive ganze Zahlen,
    deren Summe n ergibt, wobei die Reihenfolge keine Rolle spielt.

    Beispiel: p(4) = 5, denn:
        4 = 3+1 = 2+2 = 2+1+1 = 1+1+1+1

    Erzeugende Funktion (Euler):
        Σ p(n) x^n = ∏_{k≥1} 1/(1-x^k)

    @author Kurt Ingwer
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def partition_count(n: int) -> int:
        """
        Berechnet p(n) via dynamischer Programmierung.

        Rekurrenz via Euler-Pentagonalzahltheorem:
            p(n) = Σ_{k≠0} (-1)^{k+1} p(n - k(3k-1)/2)
        Alternativ: einfache DP-Tabelle mit Münz-Wechsel-Ansatz.

        Hier wird der Münz-Wechsel-Ansatz verwendet:
            dp[i] += dp[i - part]  für part = 1, 2, ..., n

        @param n: Nicht-negative ganze Zahl
        @return: Anzahl der Partitionen p(n)
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")
        if n == 0:
            return 1  # Leere Partition ist die einzige Partition von 0

        # dp[i] = Anzahl der Partitionen von i
        dp = [0] * (n + 1)
        dp[0] = 1  # Basis: eine Partition für 0 (leere Partition)

        # Füge nacheinander jeden Summanden 1, 2, ..., n hinzu
        for part in range(1, n + 1):
            for i in range(part, n + 1):
                dp[i] += dp[i - part]  # dp[i] erhöht sich um Partitionen ohne neuen Summanden

        return dp[n]

    @staticmethod
    def generate_partitions(n: int) -> list[tuple[int, ...]]:
        """
        Erzeugt alle Partitionen von n in nicht-aufsteigender Reihenfolge.

        Jede Partition wird als sortiertes Tupel (absteigend) dargestellt.
        Die Erzeugung erfolgt rekursiv: Wähle den größten Summanden ≤ n,
        dann partitioniere den Rest mit Summanden ≤ letztem Summanden.

        @param n: Nicht-negative ganze Zahl
        @return: Liste aller Partitionen als Tupel
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")

        def _gen(remaining: int, max_part: int) -> Generator:
            """Rekursiver Hilfsgenerator."""
            if remaining == 0:
                yield ()
                return
            # Wähle den aktuellen Summanden von max_part absteigend
            for part in range(min(max_part, remaining), 0, -1):
                for rest in _gen(remaining - part, part):
                    yield (part,) + rest

        return list(_gen(n, n))

    @staticmethod
    def partition_with_parts(n: int, k: int) -> list[tuple[int, ...]]:
        """
        Erzeugt alle Partitionen von n in genau k Teile.

        Eine Partition in genau k Teile bedeutet: n = a_1 + ... + a_k
        mit a_1 ≥ a_2 ≥ ... ≥ a_k ≥ 1.

        Zusammenhang: p(n, k Teile) = p(n - k, max. k Teile)
        (verschiebe jeden Summanden um 1 nach unten)

        @param n: Nicht-negative ganze Zahl
        @param k: Genau k Summanden
        @return: Liste aller Partitionen in genau k Teile
        @raises ValueError: Wenn n < 0 oder k < 0
        @lastModified 2026-03-10
        """
        if n < 0 or k < 0:
            raise ValueError("n und k müssen nicht-negativ sein")

        if k == 0:
            return [()] if n == 0 else []

        def _gen(remaining: int, max_part: int, parts_left: int) -> Generator:
            """Rekursiver Generator für Partitionen in genau k Teile."""
            if parts_left == 0:
                if remaining == 0:
                    yield ()
                return
            # Minimaler Summand: ceil(remaining / parts_left) für Gültigkeit
            min_part = max(1, -(-remaining // parts_left))  # Aufrunden ohne math.ceil
            for part in range(min(max_part, remaining - parts_left + 1), min_part - 1, -1):
                for rest in _gen(remaining - part, part, parts_left - 1):
                    yield (part,) + rest

        return list(_gen(n, n, k))

    @staticmethod
    def partition_into_distinct_parts(n: int) -> list[tuple[int, ...]]:
        """
        Erzeugt alle Partitionen von n in paarweise verschiedene Teile.

        Euler's Theorem (Symmetrie):
            Anzahl der Partitionen in verschiedene Teile
            = Anzahl der Partitionen in ungerade Teile

        Beispiel: n=6: {6, 5+1, 4+2, 3+2+1}

        @param n: Nicht-negative ganze Zahl
        @return: Liste der Partitionen in paarweise verschiedene Summanden
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")

        def _gen(remaining: int, max_part: int) -> Generator:
            """Rekursiver Generator für Partitionen in verschiedene Teile."""
            if remaining == 0:
                yield ()
                return
            # Jeder Summand ist streng kleiner als der vorherige (max_part)
            for part in range(min(max_part, remaining), 0, -1):
                for rest in _gen(remaining - part, part - 1):
                    yield (part,) + rest

        return list(_gen(n, n))

    @staticmethod
    def euler_pentagonal_theorem(n: int) -> int:
        """
        Berechnet p(n) via Euler's Pentagonalzahltheorem.

        Das Pentagonalzahltheorem (Euler, 1748):
            ∏_{k≥1}(1-x^k) = Σ_{k=-∞}^{∞} (-1)^k x^{k(3k-1)/2}

        Die verallgemeinerten Pentagonalzahlen sind:
            ω(k) = k(3k-1)/2 für k = 0, 1, -1, 2, -2, 3, -3, ...

        Die Rekurrenz für p(n):
            p(n) = p(n-1) + p(n-2) - p(n-5) - p(n-7) + p(n-12) + p(n-15) - ...

        Das Vorzeichen folgt dem Schema: +, +, -, -, +, +, -, -, ...

        @param n: Nicht-negative ganze Zahl
        @return: p(n) via Pentagonalzahl-Rekurrenz
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")
        if n == 0:
            return 1

        # Tabelle aufbauen (Bottom-Up via Pentagonalzahl-Rekurrenz)
        p = [0] * (n + 1)
        p[0] = 1

        for m in range(1, n + 1):
            total = 0
            k = 1
            # Iteration über generalisierte Pentagonalzahlen ω(k) und ω(-k)
            while True:
                # Generalisierte Pentagonalzahl: ω(k) = k(3k-1)/2
                pent_pos = k * (3 * k - 1) // 2  # k > 0
                pent_neg = k * (3 * k + 1) // 2  # -k (da (-k)(3(-k)-1)/2 = k(3k+1)/2)

                if pent_pos > m:
                    break  # Alle relevanten Pentagonalzahlen verarbeitet

                # Vorzeichen: k=1 → +, k=2 → -, k=3 → +, ...
                sign = (-1) ** (k + 1)

                total += sign * p[m - pent_pos]
                if pent_neg <= m:
                    total += sign * p[m - pent_neg]

                k += 1

            p[m] = total

        return p[n]

    @staticmethod
    def ramanujan_partition_congruences(p_vals: list[int]) -> dict:
        """
        Überprüft Ramanujans Kongruenzen für die Partitionsfunktion.

        Ramanujan (1919) entdeckte drei bemerkenswerte Kongruenzen:
            p(5k + 4)  ≡ 0 (mod 5)  für alle k ≥ 0
            p(7k + 5)  ≡ 0 (mod 7)  für alle k ≥ 0
            p(11k + 6) ≡ 0 (mod 11) für alle k ≥ 0

        Diese wurden später von Atkin und Swinnerton-Dyer vollständig bewiesen.

        @param p_vals: Liste von Werten k, für die p(5k+4), p(7k+5), p(11k+6) berechnet werden
        @return: Dict mit Verifikationsergebnissen für alle drei Kongruenzen
        @lastModified 2026-03-10
        """
        results = {
            "mod5": [],   # p(5k+4) mod 5 = 0
            "mod7": [],   # p(7k+5) mod 7 = 0
            "mod11": [],  # p(11k+6) mod 11 = 0
        }

        pf = PartitionFunction()

        for k in p_vals:
            # Kongruenz mod 5: p(5k+4) ≡ 0 (mod 5)
            n5 = 5 * k + 4
            val5 = pf.partition_count(n5)
            results["mod5"].append({
                "k": k,
                "n": n5,
                "p_n": val5,
                "remainder": val5 % 5,
                "holds": (val5 % 5 == 0)
            })

            # Kongruenz mod 7: p(7k+5) ≡ 0 (mod 7)
            n7 = 7 * k + 5
            val7 = pf.partition_count(n7)
            results["mod7"].append({
                "k": k,
                "n": n7,
                "p_n": val7,
                "remainder": val7 % 7,
                "holds": (val7 % 7 == 0)
            })

            # Kongruenz mod 11: p(11k+6) ≡ 0 (mod 11)
            n11 = 11 * k + 6
            val11 = pf.partition_count(n11)
            results["mod11"].append({
                "k": k,
                "n": n11,
                "p_n": val11,
                "remainder": val11 % 11,
                "holds": (val11 % 11 == 0)
            })

        return results


# ===========================================================================
# KLASSE 2: HARDY-RAMANUJAN-KREISMETHODE
# ===========================================================================

class HardyRamanujanCircleMethod:
    """
    Asymptotische Formeln und die Kreismethode für die Partitionsfunktion.

    Die Kreismethode (Hardy & Ramanujan, 1918; perfektioniert durch Rademacher, 1937)
    ermöglicht die exakte Berechnung von p(n) durch eine konvergente Reihe.

    Grundidee:
        p(n) = (1/2πi) ∮ F(x)/x^{n+1} dx
    wobei F(x) = ∏ 1/(1-x^k) die erzeugende Funktion ist.
    Der Integrationsweg wird in Farey-Bögen aufgeteilt.

    @author Kurt Ingwer
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def asymptotic_partition(n: int) -> float:
        """
        Hardy-Ramanujan-Asymptotik für p(n).

        Die Hauptformel (Hardy & Ramanujan, 1918):
            p(n) ~ (1 / (4n√3)) · exp(π · √(2n/3))

        Diese Formel ist für große n sehr genau (relativer Fehler → 0).
        Für n=200: p(200) = 3.972.999.029.388 (exakt),
                   Asymptotik ≈ gleiche Größenordnung.

        @param n: Positive ganze Zahl
        @return: Asymptotische Näherung für p(n)
        @raises ValueError: Wenn n ≤ 0
        @lastModified 2026-03-10
        """
        if n <= 0:
            raise ValueError(f"n muss positiv sein, erhalten: {n}")

        # Hauptterm der Hardy-Ramanujan-Asymptotik
        exponent = math.pi * math.sqrt(2 * n / 3)
        prefactor = 1.0 / (4 * n * math.sqrt(3))
        return prefactor * math.exp(exponent)

    @staticmethod
    def rademacher_formula(n: int, k_max: int = 10) -> float:
        """
        Rademacher-Formel für p(n) – exakt konvergente Summe.

        Rademacher (1937) bewies die exakte Formel:
            p(n) = (1/π√2) Σ_{k=1}^{∞} A_k(n) · √k · d/dn [sinh(π/(k) √(2/3 (n-1/24)))]
                                                               / √(n - 1/24)

        Vereinfacht mit dem Kloosterman-Summen-Term A_k(n):
            A_k(n) = Σ_{h mod k, gcd(h,k)=1} ω(h,k) · exp(-2πi·hn/k)

        Hier: Näherung mit den ersten k_max Termen (konvergiert schnell).
        Für praktische Berechnungen: n muss klein genug sein, sonst overflow.

        @param n: Positive ganze Zahl
        @param k_max: Anzahl der Summanden (Standard: 10 reicht für gute Näherung)
        @return: Näherung für p(n) via Rademacher-Reihe
        @raises ValueError: Wenn n < 1
        @lastModified 2026-03-10
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        import cmath

        def _kloosterman_sum(h: int, k: int, n: int) -> complex:
            """
            Berechnet den Dedekind-Summen-gewichteten Exponentialterm.
            Vereinfachte Version: ω(h,k) wird durch Dedekind-Summen-Näherung ersetzt.

            Für die Rademacher-Formel gilt:
                ω(h,k) = exp(πi · s(h,k))
            wobei s(h,k) die Dedekind-Summe ist.
            """
            # Dedekind-Summe s(h,k)
            def dedekind_sum(h: int, k: int) -> float:
                """s(h,k) = (1/4k) Σ_{j=1}^{k-1} cot(πj/k)·cot(πhj/k)"""
                if k <= 1:
                    return 0.0
                total = 0.0
                for j in range(1, k):
                    arg1 = math.pi * j / k
                    arg2 = math.pi * h * j / k
                    # Vermeidung von cot(0) = ∞ (j=0 ausgeschlossen)
                    s1 = math.cos(arg1) / math.sin(arg1) if abs(math.sin(arg1)) > 1e-15 else 0.0
                    s2 = math.cos(arg2) / math.sin(arg2) if abs(math.sin(arg2)) > 1e-15 else 0.0
                    total += s1 * s2
                return total / (4 * k)

            ds = dedekind_sum(h % k, k)
            omega = cmath.exp(1j * math.pi * ds)
            phase = cmath.exp(-2j * math.pi * h * n / k)
            return omega * phase

        total = 0.0
        # Äußere Summe über k = 1, 2, ..., k_max
        for k in range(1, k_max + 1):
            # A_k(n) = Σ_{h: gcd(h,k)=1} Kloosterman-Term
            A_k = sum(
                _kloosterman_sum(h, k, n)
                for h in range(1, k + 1)
                if _gcd(h, k) == 1
            )

            # Argument des sinh: π/k · √(2/3 · (n - 1/24))
            arg = math.pi / k * math.sqrt(2.0 / 3.0 * (n - 1.0 / 24.0))

            if arg > 0:
                # Ableitung von sinh(arg) / √(n-1/24): ~ cosh(arg) · (π/(k·√(2/3))) / √(n-1/24)
                d_sinh = math.cosh(arg) * math.pi / (k * math.sqrt(3.0 / 2.0)) / math.sqrt(n - 1.0 / 24.0)
            else:
                d_sinh = 0.0

            # Summand: A_k(n) · √k · d_sinh
            summand = A_k * math.sqrt(k) * d_sinh
            total += summand.real  # Nur Realteil (Ergebnis ist reell)

        # Normierung: p(n) = (1/π√2) · total
        return total / (math.pi * math.sqrt(2))

    @staticmethod
    def farey_arc_contributions(n: int, k_max: int = 5) -> list[dict]:
        """
        Berechnet die Beiträge der Farey-Bögen zur Kreismethode.

        In der Kreismethode wird der Integrationskreis in Farey-Bögen
        um die rationalen Punkte h/k (mit gcd(h,k)=1) aufgeteilt.

        Jeder Bogen um h/k trägt bei mit:
            Beitrag(h/k) ~ C(k,n) · exp(2πn · h/k · i) / k

        @param n: Positive ganze Zahl
        @param k_max: Maximale Nenner k
        @return: Liste der Beitrags-Dictionaries mit h, k und Amplitude
        @raises ValueError: Wenn n < 1
        @lastModified 2026-03-10
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

        contributions = []

        for k in range(1, k_max + 1):
            for h in range(1, k + 1):
                if _gcd(h, k) != 1:
                    continue  # Nur teilerfremde Paare (h, k)

                # Amplitude des Beitrags: exp(π√(2n/3) / k) / k² (Hauptterm)
                try:
                    amplitude = math.exp(math.pi * math.sqrt(2 * n / 3) / k) / (k * k)
                except OverflowError:
                    amplitude = float('inf')

                # Phase: exp(2πi · h · n / k)
                phase_angle = 2 * math.pi * h * n / k

                contributions.append({
                    "h": h,
                    "k": k,
                    "h_over_k": h / k,       # rationaler Punkt auf dem Einheitskreis
                    "amplitude": amplitude,   # Größe des Beitrags
                    "phase_angle": phase_angle % (2 * math.pi),  # Phase in [0, 2π)
                })

        # Sortiere nach absteigender Amplitude (wichtigste Beiträge zuerst)
        contributions.sort(key=lambda x: -x["amplitude"])
        return contributions


# ===========================================================================
# KLASSE 3: WARING-PROBLEM
# ===========================================================================

class WaringProblem:
    """
    Das Waring-Problem: Darstellung natürlicher Zahlen als Summe k-ter Potenzen.

    Edward Waring (1770) behauptete ohne Beweis:
        - Jede natürliche Zahl ist Summe von 4 Quadraten (g(2) = 4)
        - Jede natürliche Zahl ist Summe von 9 Kuben   (g(3) = 9)
        - Jede natürliche Zahl ist Summe von 19 4. Potenzen (g(4) = 19)

    Hilbert (1909) bewies die Existenz von g(k) für alle k.

    Definitionen:
        g(k) = min. Anzahl k-ter Potenzen, die jede natürliche Zahl darstellen
        G(k) = min. Anzahl k-ter Potenzen für alle hinreichend großen n

    @author Kurt Ingwer
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    # Bekannte exakte Werte von g(k) für k = 1, ..., 10
    # Quellen: Hardy & Wright, "An Introduction to the Theory of Numbers"
    _G_K_VALUES = {
        1: 1,   # g(1) = 1 (trivial: n = n)
        2: 4,   # g(2) = 4 (Lagrange, 1770)
        3: 9,   # g(3) = 9 (bewiesen 1909 durch Wieferich & Kempner)
        4: 19,  # g(4) = 19 (bewiesen 1986 durch Balasubramanian et al.)
        5: 37,  # g(5) = 37 (bewiesen 1964 durch Chen)
        6: 73,  # g(6) = 73 (bekannt)
        7: 143, # g(7) = 143 (bekannt)
        8: 279, # g(8) = 279 (bekannt)
        9: 548, # g(9) = 548 (bekannt)
        10: 1079  # g(10) = 1079 (bekannt)
    }

    # Bekannte Werte von G(k) (asymptotisch optimale Anzahl)
    _BIG_G_VALUES = {
        1: 1,   # G(1) = 1
        2: 4,   # G(2) = 4 (Lagrange)
        3: 7,   # G(3) = 7 (Linnik 1943: ≤ 7; exakt unbekannt, obere Schranke)
        4: 15,  # G(4) = 15 (Davenport, 1939)
        5: 21,  # G(5) ≤ 21 (Vaughan & Wooley)
        6: 31,  # G(6) ≤ 31
        7: 45,  # G(7) ≤ 45
        8: 62,  # G(8) ≤ 62
        9: 82,  # G(9) ≤ 82
        10: 107 # G(10) ≤ 107
    }

    @staticmethod
    def g_k(k: int) -> Optional[int]:
        """
        Gibt g(k) zurück – die minimale Anzahl k-ter Potenzen für alle n.

        g(k) ist exakt bekannt für k ≤ 10.
        Formel (für k ≥ 2, wenn die Vermutung gilt):
            g(k) = 2^k + floor((3/2)^k) - 2

        @param k: Exponent k ≥ 1
        @return: g(k) oder None wenn k > 10
        @raises ValueError: Wenn k < 1
        @lastModified 2026-03-10
        """
        if k < 1:
            raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")
        return WaringProblem._G_K_VALUES.get(k)

    @staticmethod
    def G_k(k: int) -> Optional[int]:
        """
        Gibt G(k) zurück – asymptotisch optimale Anzahl k-ter Potenzen.

        G(k) ist die minimale Anzahl r, sodass jede hinreichend große
        natürliche Zahl als Summe von r k-ten Potenzen darstellbar ist.

        Bekannte exakte Werte: G(1)=1, G(2)=4, G(4)=15.

        @param k: Exponent k ≥ 1
        @return: G(k) (bekannte Werte) oder None für k > 10
        @raises ValueError: Wenn k < 1
        @lastModified 2026-03-10
        """
        if k < 1:
            raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")
        return WaringProblem._BIG_G_VALUES.get(k)

    @staticmethod
    def represent_as_squares(n: int, k: int = 4) -> Optional[tuple[int, ...]]:
        """
        Stellt n als Summe von k Quadraten dar (Lagrange-Vier-Quadrate-Satz).

        Lagrange (1770): Jede natürliche Zahl ist Summe von ≤ 4 Quadraten.

        Algorithmus: Brute-Force mit Pruning für k ≤ 4 Quadrate.
        Für k=4: n = a² + b² + c² + d² immer lösbar (außer ggf. Ausnahmen bei k<4).

        @param n: Nicht-negative ganze Zahl
        @param k: Anzahl der Quadrate (Standard: 4)
        @return: Tupel (a₁, ..., aₖ) mit n = Σaᵢ² oder None wenn nicht darstellbar
        @raises ValueError: Wenn n < 0 oder k < 1
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")
        if k < 1:
            raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")

        def _search(remaining: int, parts_left: int, max_val: int) -> Optional[tuple]:
            """Rekursive Suche nach Quadratsummen-Darstellung."""
            if parts_left == 0:
                return () if remaining == 0 else None
            if remaining == 0:
                return (0,) * parts_left  # Rest mit Nullen auffüllen

            # Starte vom größten möglichen Quadrat
            for a in range(int(math.isqrt(min(remaining, max_val * max_val))), -1, -1):
                result = _search(remaining - a * a, parts_left - 1, a)
                if result is not None:
                    return (a,) + result
            return None

        result = _search(n, k, int(math.isqrt(n)))
        return result

    @staticmethod
    def represent_as_cubes(n: int, max_search: int = 1000) -> Optional[tuple[int, ...]]:
        """
        Stellt n als Summe von Kuben dar (Waring-Problem für k=3).

        g(3) = 9: Jede natürliche Zahl ist Summe von ≤ 9 Kuben.
        G(3) = 7: Für große n genügen 7 Kuben (Linnik 1943).

        Dieser Algorithmus sucht eine Darstellung mit möglichst wenigen Kuben
        durch Greedy-Ansatz + Backtracking.

        @param n: Nicht-negative ganze Zahl
        @param max_search: Maximale Tiefe der Suche (Standard: 1000)
        @return: Tupel von Kuben-Basen oder None wenn nicht gefunden
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")
        if n == 0:
            return (0,)

        def _is_cube(m: int) -> bool:
            """Prüft ob m ein perfekter Kubus ist."""
            if m < 0:
                return False
            c = round(m ** (1 / 3))
            return any((c + d) ** 3 == m for d in [-1, 0, 1])

        def _cube_root(m: int) -> int:
            """Ganzzahlige Kubikwurzel."""
            c = round(m ** (1 / 3))
            for d in [-1, 0, 1]:
                if (c + d) ** 3 == m:
                    return c + d
            return c

        # Greedy: subtrahere den größten Kubus, iteriere
        def _greedy(remaining: int, max_cubes: int) -> Optional[list]:
            if remaining == 0:
                return []
            if max_cubes == 0:
                return None
            # Größten Kubus ≤ remaining finden
            base = int(remaining ** (1 / 3))
            while (base + 1) ** 3 <= remaining:
                base += 1
            # Versuche verschiedene Basen (Backtracking)
            for b in range(base, 0, -1):
                rest = remaining - b ** 3
                result = _greedy(rest, max_cubes - 1)
                if result is not None:
                    return [b] + result
            return None

        # Versuche mit zunehmend mehr Kuben (1, 2, ..., 9)
        for num_cubes in range(1, min(10, max_search)):
            result = _greedy(n, num_cubes)
            if result is not None:
                return tuple(sorted(result, reverse=True))

        return None

    @staticmethod
    def four_squares_theorem(n: int) -> tuple[int, int, int, int]:
        """
        Findet eine Darstellung n = a² + b² + c² + d² (Lagrange, 1770).

        Lagrange's Vier-Quadrate-Satz: Jede nicht-negative ganze Zahl
        kann als Summe von vier Quadraten (inkl. 0²) geschrieben werden.

        Beweis-Idee: Euler's Vier-Quadrate-Identität erlaubt Multiplikation,
        und für Primzahlen existiert eine konstruktive Darstellung.

        Algorithmus hier: Direkte Suche (für moderate n effizient).

        @param n: Nicht-negative ganze Zahl
        @return: Tupel (a, b, c, d) mit n = a² + b² + c² + d², a ≥ b ≥ c ≥ d ≥ 0
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")

        # Nutze represent_as_squares mit k=4
        result = WaringProblem.represent_as_squares(n, k=4)
        if result is None:
            # Sollte nie passieren (Lagrange's Satz!)
            raise RuntimeError(f"Fehler: Vier-Quadrate-Darstellung für {n} nicht gefunden")
        # Sortiere absteigend: a ≥ b ≥ c ≥ d
        return tuple(sorted(result, reverse=True))

    @staticmethod
    def three_squares_theorem_check(n: int) -> bool:
        """
        Prüft ob n als Summe von drei Quadraten darstellbar ist (Legendre, 1798).

        Legendre-Dreiersatz: n = a² + b² + c² hat eine Lösung genau dann wenn
        n NICHT die Form 4^a(8b + 7) hat.

        Äquivalente Aussage: n ist NICHT als Summe dreier Quadrate darstellbar
        genau dann wenn n = 4^a · (8b + 7) für irgendwelche a, b ≥ 0.

        @param n: Nicht-negative ganze Zahl
        @return: True wenn n als Summe von 3 Quadraten darstellbar ist
        @raises ValueError: Wenn n < 0
        @lastModified 2026-03-10
        """
        if n < 0:
            raise ValueError(f"n muss nicht-negativ sein, erhalten: {n}")
        if n == 0:
            return True  # 0 = 0² + 0² + 0²

        m = n
        # Entferne alle Faktoren 4: m = n / 4^a
        while m % 4 == 0:
            m //= 4

        # n ist NICHT darstellbar wenn m ≡ 7 (mod 8)
        return (m % 8) != 7


# ===========================================================================
# KLASSE 4: GOLDBACH-ANALYSE
# ===========================================================================

class GoldbachAnalysis:
    """
    Analyse der Goldbach-Vermutung und verwandter additiver Probleme.

    Goldbach-Vermutung (1742): Jede gerade Zahl > 2 ist Summe zweier Primzahlen.
    Bis heute unbewiesen für alle n, aber für n ≤ 4·10^18 verifiziert.

    Schwache Goldbach-Vermutung (Vinogradov 1937, bewiesen von Helfgott 2013):
        Jede ungerade Zahl > 5 ist Summe von drei Primzahlen.

    Chen-Theorem (1966):
        Jede hinreichend große gerade Zahl = p + m, wobei m ≤ 2 Primfaktoren hat.

    @author Kurt Ingwer
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def goldbach_decompositions(n: int) -> list[tuple[int, int]]:
        """
        Findet alle Goldbach-Zerlegungen von n in zwei Primzahlen.

        Eine Goldbach-Zerlegung ist ein Paar (p, q) mit p ≤ q, p+q = n,
        wobei p und q beide prim sind.

        Goldbach-Vermutung: Diese Menge ist für jede gerade n > 2 nicht leer.

        @param n: Gerade ganze Zahl > 2
        @return: Liste aller Paare (p, q) mit p ≤ q, p+q = n, beide prim
        @raises ValueError: Wenn n ≤ 2 oder n ungerade
        @lastModified 2026-03-10
        """
        if n <= 2:
            raise ValueError(f"n muss > 2 sein, erhalten: {n}")
        if n % 2 != 0:
            raise ValueError(f"n muss gerade sein für Goldbach, erhalten: {n}")

        # Sieb bis n für effiziente Primzahlabfrage
        primes_set = set(_sieve_primes(n))

        decompositions = []
        # Suche p von 2 bis n//2 (damit p ≤ q = n-p)
        for p in range(2, n // 2 + 1):
            if p in primes_set and (n - p) in primes_set:
                decompositions.append((p, n - p))

        return decompositions

    @staticmethod
    def goldbach_count(n: int) -> int:
        """
        Zählt die Anzahl der Goldbach-Zerlegungen von n.

        Die Goldbach-Funktion r(n) = |{(p,q): p≤q, p+q=n, p,q prim}|
        ist eng verwandt mit dem Hardy-Littlewood-Hauptterm:
            r(n) ~ C₂ · ∏_{p|n, p>2} (p-1)/(p-2) · n/(ln n)²

        @param n: Gerade ganze Zahl > 2
        @return: Anzahl der Goldbach-Zerlegungen
        @raises ValueError: Wenn n ≤ 2 oder n ungerade
        @lastModified 2026-03-10
        """
        return len(GoldbachAnalysis.goldbach_decompositions(n))

    @staticmethod
    def goldbach_comet_data(n_max: int) -> list[tuple[int, int]]:
        """
        Erzeugt Goldbach-Komet-Daten: Anzahl der Zerlegungen für gerade n bis n_max.

        Der "Goldbach-Komet" ist der Graph (n, r(n)) für gerade n,
        der eine charakteristische Kometenform aufweist (dichter bei n mit
        vielen kleinen Primfaktoren).

        @param n_max: Obergrenze (inklusiv)
        @return: Liste von (n, r(n)) für alle geraden n mit 4 ≤ n ≤ n_max
        @raises ValueError: Wenn n_max < 4
        @lastModified 2026-03-10
        """
        if n_max < 4:
            raise ValueError(f"n_max muss ≥ 4 sein, erhalten: {n_max}")

        # Sieb einmal für alle n bis n_max
        primes_set = set(_sieve_primes(n_max))
        data = []

        for n in range(4, n_max + 1, 2):  # Nur gerade Zahlen ≥ 4
            count = sum(
                1 for p in range(2, n // 2 + 1)
                if p in primes_set and (n - p) in primes_set
            )
            data.append((n, count))

        return data

    @staticmethod
    def ternary_goldbach_check(n: int) -> Optional[tuple[int, int, int]]:
        """
        Prüft die schwache Goldbach-Vermutung: n = p + q + r (drei Primzahlen).

        Helfgott (2013) bewies: Jede ungerade Zahl > 5 ist Summe von 3 Primzahlen.

        Algorithmus: Für jede Primzahl p₁ ≤ n/3, prüfe ob n-p₁ als Goldbach-Zerlegung
        darstellbar ist.

        @param n: Ungerade ganze Zahl > 5
        @return: Tupel (p, q, r) oder None wenn nicht gefunden
        @raises ValueError: Wenn n ≤ 5 oder n gerade
        @lastModified 2026-03-10
        """
        if n <= 5:
            raise ValueError(f"n muss > 5 sein, erhalten: {n}")
        if n % 2 == 0:
            raise ValueError(f"n muss ungerade sein für ternäre Goldbach, erhalten: {n}")

        primes = _sieve_primes(n)
        primes_set = set(primes)

        # Für jede Primzahl p₁ prüfe ob n - p₁ als Summe zweier Primzahlen darstellbar
        for p1 in primes:
            if p1 > n // 3:
                break  # p1 ≤ p2 ≤ p3 impliziert p1 ≤ n/3
            remainder = n - p1
            # Prüfe ob remainder = p2 + p3
            for p2 in primes:
                if p2 > remainder // 2:
                    break
                p3 = remainder - p2
                if p3 in primes_set:
                    return (p1, p2, p3)  # Gefunden!

        return None  # Sollte für n > 5 ungerade nie eintreten (Helfgott)

    @staticmethod
    def chen_theorem_check(n: int) -> Optional[tuple[int, int]]:
        """
        Überprüft Chen's Theorem: n = p + m, wobei m höchstens 2 Primfaktoren hat.

        Chen Jingrun (1966) bewies: Jede hinreichend große gerade Zahl
        lässt sich als Summe einer Primzahl und einer Semiprimzahl (oder Primzahl)
        schreiben.

        Semiprimzahl: m = p·q (Produkt genau zweier Primzahlen, evtl. gleich).

        @param n: Gerade ganze Zahl > 2
        @return: Tupel (p, m) mit p prim, m Primzahl oder Semiprimzahl, p+m=n
        @raises ValueError: Wenn n ≤ 2 oder n ungerade
        @lastModified 2026-03-10
        """
        if n <= 2:
            raise ValueError(f"n muss > 2 sein, erhalten: {n}")
        if n % 2 != 0:
            raise ValueError(f"n muss gerade sein, erhalten: {n}")

        def _count_prime_factors(m: int) -> int:
            """Zählt die Primfaktoren mit Vielfachheit (Omega-Funktion Ω(m))."""
            if m < 2:
                return 0
            count = 0
            d = 2
            while d * d <= m:
                while m % d == 0:
                    count += 1
                    m //= d
                d += 1
            if m > 1:
                count += 1
            return count

        primes_set = set(_sieve_primes(n))

        for p in sorted(primes_set):
            if p >= n:
                break
            m = n - p
            if m >= 2 and _count_prime_factors(m) <= 2:
                return (p, m)  # p + m = n, m hat ≤ 2 Primfaktoren

        return None


# ===========================================================================
# KLASSE 5: ARITHMETISCHE PROGRESSIONEN
# ===========================================================================

class ArithmeticProgressions:
    """
    Arithmetische Progressionen in den Primzahlen und verwandten Strukturen.

    Wichtige Sätze:
    - Dirichlet (1837): Es gibt unendlich viele Primzahlen in jeder AP a, a+d, a+2d, ...
      mit gcd(a,d) = 1.
    - Green-Tao (2004): Die Primzahlen enthalten arithmetische Progressionen
      beliebiger Länge.
    - Van-der-Waerden (1927): Für jede Färbung {1,...,N} mit r Farben
      gibt es eine einfarbige AP der Länge k (wenn N ≥ W(k;r)).

    @author Kurt Ingwer
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    # Bekannte Van-der-Waerden-Zahlen W(k; r) für kleine k, r
    # W(k; r): kleinste N, sodass jede r-Färbung von {1,...,N} eine einfarbige k-AP enthält
    # Format: (k, r) → W(k, r)
    _VAN_DER_WAERDEN = {
        (3, 2): 9,     # W(3;2) = 9 (bekannt exakt)
        (4, 2): 35,    # W(4;2) = 35
        (5, 2): 178,   # W(5;2) = 178 (Chvátal, 1970)
        (6, 2): 1132,  # W(6;2) ≥ 1132 (untere Schranke)
        (3, 3): 27,    # W(3;3) = 27
        (3, 4): 76,    # W(3;4) ≥ 76
        (4, 3): 293,   # W(4;3) ≥ 293
    }

    @staticmethod
    def primes_in_arithmetic_progression(a: int, d: int, n_max: int) -> list[int]:
        """
        Findet alle Primzahlen ≤ n_max in der arithmetischen Progression a, a+d, a+2d, ...

        Dirichlet-Satz (1837): Falls gcd(a, d) = 1, enthält die AP unendlich viele Primzahlen.
        Die Dichte beträgt asymptotisch 1/φ(d) aller Primzahlen.

        @param a: Startelement (a ≥ 1)
        @param d: Differenz (d ≥ 1)
        @param n_max: Obergrenze
        @return: Liste der Primzahlen in der AP, sortiert aufsteigend
        @raises ValueError: Wenn a < 1 oder d < 1
        @lastModified 2026-03-10
        """
        if a < 1:
            raise ValueError(f"a muss ≥ 1 sein, erhalten: {a}")
        if d < 1:
            raise ValueError(f"d muss ≥ 1 sein, erhalten: {d}")

        primes_set = set(_sieve_primes(n_max))

        # Generiere alle Glieder der AP und filtere Primzahlen
        result = []
        term = a
        while term <= n_max:
            if term in primes_set:
                result.append(term)
            term += d

        return result

    @staticmethod
    def van_der_waerden_w(k: int, r: int) -> Optional[int]:
        """
        Gibt die Van-der-Waerden-Zahl W(k; r) zurück (soweit bekannt).

        W(k; r) ist die kleinste natürliche Zahl N, sodass jede r-Färbung
        von {1, ..., N} eine einfarbige arithmetische Progression der Länge k
        enthält.

        Van-der-Waerden-Satz (1927): W(k; r) existiert für alle k, r ≥ 1.

        @param k: AP-Länge k ≥ 2
        @param r: Anzahl der Farben r ≥ 2
        @return: W(k; r) wenn bekannt, sonst None
        @raises ValueError: Wenn k < 2 oder r < 2
        @lastModified 2026-03-10
        """
        if k < 2:
            raise ValueError(f"k muss ≥ 2 sein, erhalten: {k}")
        if r < 2:
            raise ValueError(f"r muss ≥ 2 sein, erhalten: {r}")
        return ArithmeticProgressions._VAN_DER_WAERDEN.get((k, r))

    @staticmethod
    def arithmetic_progression_in_primes(length: int, start_search: int = 2) -> Optional[tuple]:
        """
        Findet eine arithmetische Progression der Länge `length` in den Primzahlen.

        Green-Tao-Theorem (2004): Die Menge der Primzahlen enthält
        arithmetische Progressionen beliebiger endlicher Länge.

        Algorithmus: Suche durch alle Paare (p, d) bis zu einem vernünftigen Limit.

        @param length: Gewünschte AP-Länge (≥ 2)
        @param start_search: Startpunkt der Suche (Standard: 2)
        @return: Tupel (p, d) sodass p, p+d, ..., p+(length-1)d alle prim sind, oder None
        @raises ValueError: Wenn length < 2
        @lastModified 2026-03-10
        """
        if length < 2:
            raise ValueError(f"length muss ≥ 2 sein, erhalten: {length}")

        # Suchgrenze: heuristisch angepasst an die AP-Länge
        search_limit = max(1000, 500 * length * length)
        primes = _sieve_primes(search_limit)
        primes_set = set(primes)

        # Durchsuche alle möglichen Startprimzahlen und Differenzen
        for p in primes:
            if p < start_search:
                continue
            # Differenz d: maximal so groß, dass p + (length-1)*d ≤ search_limit
            max_d = (search_limit - p) // (length - 1)
            for d in range(2, max_d + 1, 2):  # d immer gerade für length > 2 (außer p=2)
                # Prüfe ob p, p+d, p+2d, ..., p+(length-1)d alle prim sind
                ap = [p + i * d for i in range(length)]
                if all(q in primes_set for q in ap):
                    return tuple(ap)

        return None

    @staticmethod
    def dirichlet_density(a: int, d: int) -> float:
        """
        Berechnet die Dirichlet-Dichte der Primzahlen in der AP a mod d.

        Dirichlet-Satz (1837) mit Dichteaussage:
            Für gcd(a, d) = 1 gilt: d(P_{a,d}) = 1/φ(d)

        Das bedeutet: Der Anteil der Primzahlen in der AP a, a+d, a+2d, ...
        an allen Primzahlen beträgt asymptotisch 1/φ(d).

        @param a: Startglied
        @param d: Differenz
        @return: Asymptotische Dichte 1/φ(d) wenn gcd(a,d)=1, sonst 0
        @raises ValueError: Wenn d < 1
        @lastModified 2026-03-10
        """
        if d < 1:
            raise ValueError(f"d muss ≥ 1 sein, erhalten: {d}")

        # Prüfe ob gcd(a, d) = 1 (Voraussetzung für Dirichlet)
        if _gcd(a % d, d) != 1:
            return 0.0  # Nur endlich viele Primzahlen wenn gcd ≠ 1

        # Euler'sche Phi-Funktion φ(d) berechnen
        def _euler_phi(n: int) -> int:
            """φ(n) = n · ∏_{p|n} (1 - 1/p)"""
            result = n
            temp = n
            p = 2
            while p * p <= temp:
                if temp % p == 0:
                    while temp % p == 0:
                        temp //= p
                    result -= result // p  # result *= (1 - 1/p)
                p += 1
            if temp > 1:
                result -= result // temp
            return result

        phi_d = _euler_phi(d)
        return 1.0 / phi_d

    @staticmethod
    def green_tao_demo() -> dict:
        """
        Demonstriert das Green-Tao-Theorem anhand konkreter APs in Primzahlen.

        Green-Tao-Theorem (Terence Tao, Ben Green, 2004):
            Die Menge der Primzahlen enthält arithmetische Progressionen
            beliebiger endlicher Länge.

        Dieser Beweis war eine der größten Errungenschaften in der Zahlentheorie
        des 21. Jahrhunderts und nutzt Techniken aus der Ergoden-Theorie.

        Bekannte lange APs in Primzahlen:
            - Länge 6: 7, 157, 307, 457, 607, 757 (d=150)
            - Länge 7: 7, 157, 307, 457, 607, 757, 907 (d=150)
            - Länge 21: erste bekannte AP begann mit p=...

        @return: Dictionary mit bekannten Beispielen und Erläuterungen
        @lastModified 2026-03-10
        """
        # Verifiziere bekannte Beispiele
        examples = []

        # Beispiel 1: AP der Länge 3: 3, 7, 11 (d=4)
        ap3 = (3, 7, 11)
        examples.append({
            "length": 3,
            "progression": ap3,
            "difference": 4,
            "all_prime": all(_is_prime_simple(p) for p in ap3)
        })

        # Beispiel 2: AP der Länge 5: 5, 11, 17, 23, 29 (d=6)
        ap5 = (5, 11, 17, 23, 29)
        examples.append({
            "length": 5,
            "progression": ap5,
            "difference": 6,
            "all_prime": all(_is_prime_simple(p) for p in ap5)
        })

        # Beispiel 3: AP der Länge 6: 7, 157, 307, 457, 607, 757 (d=150)
        ap6 = (7, 157, 307, 457, 607, 757)
        examples.append({
            "length": 6,
            "progression": ap6,
            "difference": 150,
            "all_prime": all(_is_prime_simple(p) for p in ap6)
        })

        return {
            "theorem": "Green-Tao (2004): Primzahlen enthalten APs beliebiger Länge",
            "examples": examples,
            "note": "Bewiesen via Ergoden-Theorie und Fourier-Analyse",
            "year": 2004,
            "fields_medal": 2006  # Tao erhielt die Fields-Medaille für diesen Beweis
        }


# ===========================================================================
# KLASSE 6: SUMMENMENGENTHEORIE
# ===========================================================================

class SumSetTheory:
    """
    Summenmengen und additive Kombinatorik.

    Die additive Kombinatorik untersucht die Struktur von Mengen A, B ⊆ ℤ
    bezüglich ihrer Summen- und Differenzmengen.

    Wichtige Resultate:
    - Cauchy-Davenport: |A+B| ≥ min(p, |A|+|B|-1) für Mengen mod p
    - Plünnecke-Ruzsa: |kA| / |A| ≤ (|A+B|/|A|)^k
    - Freiman's Theorem: Kleine Verdoppelungskonstante → A liegt in einem GAP

    @author Kurt Ingwer
    @since 2026-03-10
    @lastModified 2026-03-10
    """

    @staticmethod
    def sumset(A: set, B: set) -> set:
        """
        Berechnet die Summenmenge A + B = {a + b | a ∈ A, b ∈ B}.

        Die Summenmenge ist ein zentrales Objekt der additiven Kombinatorik.
        Für Intervalle gilt: |[1,n] + [1,m]| = n + m - 1.

        @param A: Erste Menge ganzer Zahlen
        @param B: Zweite Menge ganzer Zahlen
        @return: Summenmenge A + B
        @lastModified 2026-03-10
        """
        if not A or not B:
            return set()  # Leere Menge + beliebige Menge = leere Menge
        return {a + b for a in A for b in B}

    @staticmethod
    def difference_set(A: set, B: set) -> set:
        """
        Berechnet die Differenzmenge A - B = {a - b | a ∈ A, b ∈ B}.

        Die Differenzmenge erfüllt: |A - A| ≥ 2|A| - 1 (Cauchy).
        Für Mengen mit kleinem Quotient |A-A|/|A| gilt Freiman's Theorem.

        @param A: Erste Menge ganzer Zahlen (Minuend)
        @param B: Zweite Menge ganzer Zahlen (Subtrahend)
        @return: Differenzmenge A - B
        @lastModified 2026-03-10
        """
        if not A or not B:
            return set()  # Leere Differenzmenge
        return {a - b for a in A for b in B}

    @staticmethod
    def cauchy_davenport_bound(A: set, B: set, p: int) -> int:
        """
        Berechnet die Cauchy-Davenport-Schranke |A+B| ≥ min(p, |A|+|B|-1).

        Cauchy-Davenport-Theorem (1813/1935):
            Für nichtleere Mengen A, B ⊆ ℤ/pℤ (p Primzahl) gilt:
            |A + B| ≥ min(p, |A| + |B| - 1)

        Dies ist die schärfste untere Schranke für |A+B| in ℤ/pℤ.
        Gleichheit gilt genau wenn A oder B ein arithmetisches Intervall ist.

        @param A: Erste Menge (Elemente als ganze Zahlen, werden mod p interpretiert)
        @param B: Zweite Menge
        @param p: Primzahl (Modulus)
        @return: Die Cauchy-Davenport-Schranke min(p, |A|+|B|-1)
        @raises ValueError: Wenn p keine Primzahl oder Mengen leer
        @lastModified 2026-03-10
        """
        if not A or not B:
            raise ValueError("A und B müssen nicht leer sein")
        if not _is_prime_simple(p):
            raise ValueError(f"p muss eine Primzahl sein, erhalten: {p}")

        # Reduziere A und B modulo p
        A_mod = {a % p for a in A}
        B_mod = {b % p for b in B}

        # Cauchy-Davenport-Schranke
        return min(p, len(A_mod) + len(B_mod) - 1)

    @staticmethod
    def plunnecke_ruzsa_estimate(A: set, B: set) -> dict:
        """
        Berechnet die Plünnecke-Ruzsa-Ungleichung für Summenmengen.

        Plünnecke-Ruzsa-Ungleichung:
            Falls K = |A+B|/|A| (Verdoppelungskonstante bezgl. B),
            dann gilt für k, l ≥ 0:
            |kB - lB| ≤ K^{k+l} · |A|

        Spezialfall (Ruzsa-Dreiecksungleichung):
            |A - C| ≤ |A - B| · |B - C| / |B|

        @param A: Referenzmenge
        @param B: Summationsmenge
        @return: Dictionary mit Verdoppelungskonstante und Schranken
        @raises ValueError: Wenn A oder B leer
        @lastModified 2026-03-10
        """
        if not A or not B:
            raise ValueError("A und B dürfen nicht leer sein")

        # Summenmenge A + B
        sumset_AB = SumSetTheory.sumset(A, B)
        # Differenzmenge A - A (Ruzsa-Distanz)
        diffset_AA = SumSetTheory.difference_set(A, A)

        # Verdoppelungskonstante K = |A+B| / |A|
        K = len(sumset_AB) / len(A)

        # Schranken für höhere Iterated Sumsets via Plünnecke
        # |A + kB| ≤ K^k · |A| (obere Schranke)
        bounds = {}
        for k in range(1, 5):
            bounds[f"|A + {k}B|_upper"] = K ** k * len(A)

        return {
            "K_doubling": K,                      # Verdoppelungskonstante
            "|A+B|": len(sumset_AB),              # tatsächliche Summenmengen-Größe
            "|A-A|": len(diffset_AA),             # Differenzmengen-Größe
            "ruzsa_quotient": len(diffset_AA) / len(A),  # Ruzsa-Quotient
            "plunnecke_bounds": bounds,           # Plünnecke-Schranken
        }

    @staticmethod
    def freiman_theorem_demo(A: set) -> dict:
        """
        Demonstration von Freiman's Theorem für eine Menge A.

        Freiman's Theorem (1964):
            Wenn |A + A| ≤ K · |A| für ein kleines K,
            dann liegt A in einem verallgemeinerten arithmetischen Progressionen (GAP)
            mit Dimension d ≤ d(K) und Größe ≤ f(K) · |A|.

        Für K = 2: A ist (fast) ein arithmetisches Intervall.
        Für K < 3/2: A muss in einer AP der Richtigen Größe enthalten sein.

        @param A: Menge ganzer Zahlen
        @return: Dictionary mit Verdoppelungskonstante und Strukturanalyse
        @raises ValueError: Wenn A leer
        @lastModified 2026-03-10
        """
        if not A:
            raise ValueError("A darf nicht leer sein")

        # Summenmenge A + A (Iterated Sumset der Ordnung 2)
        sumset_2A = SumSetTheory.sumset(A, A)
        K = len(sumset_2A) / len(A)

        # Klassifikation nach Freiman
        if K < 2:
            structure = "A liegt in einer arithmetischen Progression (AP)"
        elif K < 3:
            structure = "A liegt in einer generalisierten AP (GAP) kleiner Dimension"
        elif K < 4:
            structure = "A hat mittlere additive Struktur"
        else:
            structure = "A hat wenig additive Struktur (Freiman's Theorem gibt großes GAP)"

        # Überprüfe ob A selbst eine AP ist
        sorted_A = sorted(A)
        is_ap = len(A) <= 1 or all(
            sorted_A[i + 1] - sorted_A[i] == sorted_A[1] - sorted_A[0]
            for i in range(len(sorted_A) - 1)
        )

        return {
            "K_doubling": K,          # Verdoppelungskonstante |A+A|/|A|
            "|A+A|": len(sumset_2A),  # Größe der iterierten Summenmenge
            "|A|": len(A),            # Größe der Ausgangsmenge
            "is_arithmetic_progression": is_ap,  # Ist A selbst eine AP?
            "structure_classification": structure,
        }

    @staticmethod
    def additive_energy(A: set) -> int:
        """
        Berechnet die additive Energie E(A) einer Menge.

        Definition:
            E(A) = |{(a, b, c, d) ∈ A⁴ : a + b = c + d}|

        Die additive Energie misst die "additive Struktur" von A.
        Es gilt: |A|² ≤ E(A) ≤ |A|³ (triviale Schranken).

        Zusammenhang mit Summenmengen (Cauchy-Schwarz):
            |A + A| ≥ |A|⁴ / E(A)

        Für Mengen mit kleinem |A+A| ist E(A) groß (Satz von Balog-Szemerédi-Gowers).

        @param A: Menge ganzer Zahlen
        @return: Additive Energie E(A)
        @raises ValueError: Wenn A leer
        @lastModified 2026-03-10
        """
        if not A:
            raise ValueError("A darf nicht leer sein")

        A_list = sorted(A)
        n = len(A_list)

        # Zähle Lösungen von a + b = c + d
        # Effizient via Häufigkeiten der Differenzen (Faltung)
        # E(A) = Σ_{s} r(s)² wobei r(s) = |{(a,b): a-b=s, a,b∈A}|

        from collections import Counter
        # Zähle Häufigkeiten der Differenzen a - b
        diff_count: Counter = Counter()
        for a in A_list:
            for b in A_list:
                diff_count[a - b] += 1

        # E(A) = Σ r(s)²
        energy = sum(count * count for count in diff_count.values())
        return energy


# ===========================================================================
# STANDALONE-FUNKTIONEN
# ===========================================================================

# Bekannte Schur-Zahlen S(k): Kleinste N, sodass jede k-Färbung von {1,...,N}
# eine monochromatische Lösung von a + b = c enthält.
_SCHUR_NUMBERS = {
    1: 2,    # S(1) = 2
    2: 5,    # S(2) = 5
    3: 14,   # S(3) = 14
    4: 45,   # S(4) = 45
    5: 161,  # S(5) = 161 (bewiesen 2017 via SAT-Solver: Heule et al.)
}


def schur_number(k: int) -> Optional[int]:
    """
    Gibt die Schur-Zahl S(k) zurück.

    Schur-Satz (1916): Für jedes k gibt es eine Zahl S(k), sodass
    jede k-Färbung von {1, ..., S(k)} eine monochromatische Lösung
    der Gleichung a + b = c enthält.

    Bekannte Werte:
        S(1) = 2, S(2) = 5, S(3) = 14, S(4) = 45, S(5) = 161

    S(5) = 161 wurde 2017 durch einen 2-Petabyte-SAT-Beweis (Heule et al.)
    bestätigt – einer der größten Computerbeweis der Geschichte.

    @param k: Anzahl der Farben k ≥ 1
    @return: S(k) wenn k ≤ 5, sonst None (unbekannt)
    @raises ValueError: Wenn k < 1
    @lastModified 2026-03-10
    """
    if k < 1:
        raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")
    return _SCHUR_NUMBERS.get(k)


def happy_ending_problem(n: int) -> int:
    """
    Berechnet die Erdős-Szekeres-Schranke für das Happy-Ending-Problem.

    Happy-Ending-Problem (Esther Klein, 1935):
        Unter je 5 Punkten in allgemeiner Lage (keine 3 kollinear)
        gibt es immer 4 Punkte, die ein konvexes Viereck bilden.

    Erdős-Szekeres-Satz (1935):
        Für je f(n) = C(2n-4, n-2) + 1 Punkte gibt es immer n Punkte
        die ein konvexes n-Eck bilden.

    Allgemein: f(n) ≤ 2^{n-2} + 1 (obere Schranke, Erdős-Szekeres)

    @param n: Größe des konvexen Polygons (n ≥ 3)
    @return: Erdős-Szekeres-Schranke: kleinste Punktzahl garantierend n-Eck
    @raises ValueError: Wenn n < 3
    @lastModified 2026-03-10
    """
    if n < 3:
        raise ValueError(f"n muss ≥ 3 sein (konvexes Polygon), erhalten: {n}")

    # Erdős-Szekeres-Schranke: f(n) = 2^{n-2} + 1
    # Für n=4: f(4) = 5 (bestätigt durch Esther Kleins Beobachtung)
    # Für n=5: f(5) = 9 (bewiesen)
    # Für n≥6: 2^{n-2}+1 (obere Schranke, exakter Wert unbekannt)
    return 2 ** (n - 2) + 1


def bertrand_postulate_verify(n: int) -> Optional[int]:
    """
    Verifiziert Bertrand's Postulat: Zwischen n und 2n liegt immer eine Primzahl.

    Bertrand's Postulat (1845, bewiesen von Tschebyschow 1852):
        Für jedes n ≥ 1 gibt es eine Primzahl p mit n < p ≤ 2n.

    Ramanujan (1919) und Erdős (1932) lieferten weitere elementare Beweise.

    @param n: Positive ganze Zahl n ≥ 1
    @return: Eine Primzahl p mit n < p ≤ 2n
    @raises ValueError: Wenn n < 1
    @lastModified 2026-03-10
    """
    if n < 1:
        raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")

    # Suche die kleinste Primzahl im offenen Intervall (n, 2n]
    for p in range(n + 1, 2 * n + 1):
        if _is_prime_simple(p):
            return p  # Gefunden (Bertrand's Postulat garantiert Existenz)

    return None  # Sollte für n ≥ 1 nie eintreten (Tschebyschow)
