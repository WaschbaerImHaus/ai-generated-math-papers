"""
@file combinatorics.py
@brief Kombinatorik und Diskrete Mathematik – Permutationsgruppen, abzählende Kombinatorik,
       erzeugende Funktionen, Ramsey-Theorie, Graphen-Kombinatorik, Kombinatorik auf Wörtern.
@author Kurt Ingwer
@lastModified 2026-03-10

Mathematische Grundlagen:
    - Permutationsgruppe S_n: alle bijektiven Abbildungen {1,...,n} → {1,...,n}
    - Stirling-Zahlen 1. Art s(n,k): Anzahl Permutationen mit k Zyklen (vorzeichenbehaftet)
    - Stirling-Zahlen 2. Art S(n,k): Anzahl Partitionen einer n-Menge in k nicht-leere Blöcke
    - Ramsey-Zahl R(s,t): kleinste n, sodass jede 2-Färbung von K_n
      einen roten K_s oder blauen K_t enthält
    - Kirchhoff-Matrix-Satz: Anzahl Spannbäume = det(L̃) wobei L̃ reduzierte Laplace-Matrix
"""

import math
import itertools
from math import gcd, factorial, comb
from functools import lru_cache
from typing import Any


# ============================================================
# Hilfsfunktionen
# ============================================================

def _lcm(a: int, b: int) -> int:
    """
    @brief Kleinstes gemeinsames Vielfaches von a und b.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    @param a Erste ganze Zahl.
    @param b Zweite ganze Zahl.
    @return kgV(a, b).
    """
    # kgV = |a·b| / ggT(a,b)
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)


# ============================================================
# 1. PermutationGroup – Permutationsgruppen
# ============================================================

class PermutationGroup:
    """
    @brief Symmetrische Gruppe S_n – alle Permutationen von {1, ..., n}.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Eine Permutation wird als Liste dargestellt: perm[i] = Bild von i (0-indiziert).
    Beispiel: perm = [2, 0, 1] bedeutet 0→2, 1→0, 2→1.

    Mathematisch: S_n = Gruppe aller Bijektionen {0,...,n-1} → {0,...,n-1}
    mit Komposition als Gruppenoperation. |S_n| = n!
    """

    def __init__(self, n: int):
        """
        @brief Initialisiert die symmetrische Gruppe S_n.
        @author Kurt Ingwer
        @lastModified 2026-03-10
        @param n Grad der symmetrischen Gruppe (Anzahl der Elemente).
        """
        # Grad der Gruppe speichern
        self.n = n
        # Identitätspermutation: jedes Element bleibt an seiner Stelle
        self.identity = list(range(n))

    @staticmethod
    def all_permutations(n: int) -> list[list[int]]:
        """
        @brief Erzeugt alle Permutationen von {0, ..., n-1}.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Mathematisch: Auflistung aller Elemente der symmetrischen Gruppe S_n.
        |S_n| = n! Elemente.

        @param n Anzahl der Elemente.
        @return Liste aller n! Permutationen als Listen.
        """
        # itertools.permutations liefert alle möglichen Anordnungen
        return [list(p) for p in itertools.permutations(range(n))]

    @staticmethod
    def permutation_sign(perm: list[int]) -> int:
        """
        @brief Berechnet das Vorzeichen (Signum) einer Permutation.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Das Vorzeichen ist +1 (gerade) oder -1 (ungerade), je nach Anzahl der
        Transpositionen in der Zerlegung der Permutation.

        Methode: Zähle Inversionen (Paare i<j mit perm[i]>perm[j]).
        sign(σ) = (-1)^{Anzahl Inversionen}

        @param perm Permutation als Liste.
        @return +1 für gerade, -1 für ungerade Permutation.
        """
        n = len(perm)
        if n == 0:
            return 1  # leere Permutation ist gerade (0 Inversionen)

        # Anzahl der Inversionen zählen
        inversions = 0
        for i in range(n):
            for j in range(i + 1, n):
                # Inversion: Element an Position i ist größer als an Position j
                if perm[i] > perm[j]:
                    inversions += 1

        # Gerade Anzahl → +1, ungerade → -1
        return 1 if inversions % 2 == 0 else -1

    @staticmethod
    def cycle_decomposition(perm: list[int]) -> list[list[int]]:
        """
        @brief Zerlegt eine Permutation in disjunkte Zyklen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Jede Permutation lässt sich eindeutig als Produkt disjunkter Zyklen
        schreiben (bis auf Reihenfolge und Zyklenrotation).

        Beispiel: [1,2,0,4,3] → [[0,1,2], [3,4]] (3-Zyklus und 2-Zyklus)

        @param perm Permutation als Liste (0-indiziert).
        @return Liste der Zyklen (Fixpunkte als 1-elementige Listen eingeschlossen).
        """
        n = len(perm)
        visited = [False] * n   # Welche Elemente wurden schon einem Zyklus zugeordnet?
        cycles = []

        for start in range(n):
            if visited[start]:
                continue  # Bereits in einem Zyklus

            # Neuen Zyklus von 'start' aus verfolgen
            cycle = []
            current = start
            while not visited[current]:
                visited[current] = True
                cycle.append(current)
                current = perm[current]  # Nächstes Element im Zyklus

            cycles.append(cycle)

        return cycles

    @staticmethod
    def order_of_permutation(perm: list[int]) -> int:
        """
        @brief Berechnet die Ordnung einer Permutation.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Ordnung von σ ∈ S_n ist die kleinste positive Zahl k mit σ^k = id.
        Sie entspricht dem kgV der Längen aller Zyklen in der Zyklenzerlegung.

        Beweis: σ^k = id ⟺ k ist Vielfaches jeder Zyklenlänge.

        @param perm Permutation als Liste.
        @return Ordnung der Permutation (positive ganze Zahl).
        """
        if len(perm) == 0:
            return 1  # leere Permutation: Ordnung 1

        # Zyklenzerlegung berechnen
        cycles = PermutationGroup.cycle_decomposition(perm)

        # Ordnung = kgV aller Zyklenlängen
        result = 1
        for cycle in cycles:
            result = _lcm(result, len(cycle))

        return result

    @staticmethod
    def compose_permutations(p1: list[int], p2: list[int]) -> list[int]:
        """
        @brief Berechnet die Komposition p1 ∘ p2 (erst p2, dann p1).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        (p1 ∘ p2)(i) = p1(p2(i))

        @param p1 Erste Permutation (äußere).
        @param p2 Zweite Permutation (innere).
        @return Komposition p1 ∘ p2 als Liste.
        @raises ValueError bei unterschiedlicher Länge.
        """
        if len(p1) != len(p2):
            raise ValueError("Permutationen müssen gleiche Länge haben.")

        # Komposition: erst p2 anwenden, dann p1
        return [p1[p2[i]] for i in range(len(p1))]

    @staticmethod
    def inverse_permutation(p: list[int]) -> list[int]:
        """
        @brief Berechnet die inverse Permutation p^{-1}.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        p^{-1} ist definiert durch p(p^{-1}(i)) = i für alle i.
        Äquivalent: inv[p[i]] = i.

        @param p Permutation als Liste.
        @return Inverse Permutation als Liste.
        """
        n = len(p)
        inv = [0] * n

        # inv[p[i]] = i: Umkehrabbildung berechnen
        for i in range(n):
            inv[p[i]] = i

        return inv

    @staticmethod
    def is_even(perm: list[int]) -> bool:
        """
        @brief Prüft, ob eine Permutation gerade ist (zur alternierenden Gruppe A_n gehört).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        A_n = {σ ∈ S_n | sign(σ) = +1} ist Untergruppe vom Index 2 in S_n.
        |A_n| = n!/2 für n ≥ 2.

        @param perm Permutation als Liste.
        @return True wenn gerade (sign = +1), False wenn ungerade.
        """
        return PermutationGroup.permutation_sign(perm) == 1

    @staticmethod
    def fixed_points(perm: list[int]) -> list[int]:
        """
        @brief Gibt alle Fixpunkte der Permutation zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein Fixpunkt i erfüllt perm[i] = i (wird nicht bewegt).
        Mathematisch: Fix(σ) = {i ∈ {0,...,n-1} | σ(i) = i}.

        @param perm Permutation als Liste.
        @return Liste der Fixpunkte.
        """
        # Elemente suchen, die auf sich selbst abgebildet werden
        return [i for i in range(len(perm)) if perm[i] == i]


# ============================================================
# 2. EnumerativeCombinatorics – Abzählende Kombinatorik
# ============================================================

class EnumerativeCombinatorics:
    """
    @brief Abzählende Kombinatorik – Multinomialkoeffizienten, Stirling-Zahlen,
           Euler-Zahlen, Partitionszahlen und weitere kombinatorische Kennzahlen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Stirling 1. Art: s(n,k) = s(n-1,k-1) - (n-1)·s(n-1,k), s(0,0)=1
    Stirling 2. Art: S(n,k) = k·S(n-1,k) + S(n-1,k-1), S(0,0)=1
    Partition: p(n) = Anzahl Wege n als Summe positiver ganzer Zahlen zu schreiben
    """

    @staticmethod
    def multinomial(n: int, ks: list[int]) -> int:
        """
        @brief Berechnet den Multinomialkoeffizienten n!/(k1! · k2! · ... · km!).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Gibt die Anzahl Möglichkeiten an, n Objekte in Gruppen der Größen k1,...,km aufzuteilen.
        Bedingung: k1 + k2 + ... + km = n.

        @param n Gesamtanzahl der Objekte.
        @param ks Liste der Gruppengrößen (müssen sich zu n addieren).
        @return Multinomialkoeffizient als ganze Zahl.
        @raises ValueError wenn sum(ks) ≠ n.
        """
        if sum(ks) != n:
            raise ValueError(f"Summe der Gruppengrößen muss n={n} ergeben, ist aber {sum(ks)}.")

        # Multinomialkoeffizient = n! / (k1! * k2! * ... * km!)
        result = factorial(n)
        for k in ks:
            result //= factorial(k)
        return result

    @staticmethod
    @lru_cache(maxsize=None)
    def stirling_first(n: int, k: int) -> int:
        """
        @brief Berechnet die vorzeichenbehaftete Stirling-Zahl 1. Art s(n,k).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        s(n,k) = Anzahl der Permutationen in S_n mit genau k Zyklen (vorzeichenbehaftet).
        Vorzeichen: s(n,k) = (-1)^{n-k} · |s(n,k)|

        Rekurrenz: s(n,k) = s(n-1,k-1) - (n-1)·s(n-1,k)
        Randbedingungen: s(0,0)=1, s(n,0)=0 für n>0, s(0,k)=0 für k>0.

        @param n Anzahl der Elemente.
        @param k Anzahl der Zyklen.
        @return Vorzeichenbehaftete Stirling-Zahl 1. Art.
        """
        # Randbedingungen
        if n == 0 and k == 0:
            return 1
        if n == 0 or k == 0:
            return 0
        if k > n:
            return 0

        # Rekurrenz: s(n,k) = s(n-1,k-1) - (n-1)*s(n-1,k)
        return (EnumerativeCombinatorics.stirling_first(n - 1, k - 1)
                - (n - 1) * EnumerativeCombinatorics.stirling_first(n - 1, k))

    @staticmethod
    @lru_cache(maxsize=None)
    def stirling_second(n: int, k: int) -> int:
        """
        @brief Berechnet die Stirling-Zahl 2. Art S(n,k).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        S(n,k) = Anzahl der Partitionen einer n-elementigen Menge in k nicht-leere Teilmengen.

        Rekurrenz: S(n,k) = k·S(n-1,k) + S(n-1,k-1)
        Randbedingungen: S(0,0)=1, S(n,0)=0 für n>0, S(0,k)=0 für k>0.

        Explizite Formel: S(n,k) = (1/k!) Σ_{j=0}^{k} (-1)^{k-j} C(k,j) j^n

        @param n Anzahl der Elemente.
        @param k Anzahl der Blöcke.
        @return Stirling-Zahl 2. Art (nicht-negativ).
        """
        # Randbedingungen
        if n == 0 and k == 0:
            return 1
        if n == 0 or k == 0:
            return 0
        if k > n:
            return 0

        # Rekurrenz: S(n,k) = k*S(n-1,k) + S(n-1,k-1)
        return (k * EnumerativeCombinatorics.stirling_second(n - 1, k)
                + EnumerativeCombinatorics.stirling_second(n - 1, k - 1))

    @staticmethod
    @lru_cache(maxsize=None)
    def euler_number_a(n: int, k: int) -> int:
        """
        @brief Berechnet die Euler-Zahlen A(n,k) – Anzahl der Permutationen von {1,...,n}
               mit genau k Anstiegen (Aufsteige: σ(i) < σ(i+1)).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Rekurrenz: A(n,k) = (k+1)·A(n-1,k) + (n-k)·A(n-1,k-1)
        Randbedingungen: A(0,0)=1, A(n,k)=0 für k<0 oder k>=n.

        Hinweis: Nicht zu verwechseln mit Euler-Zahlen aus der Analysis.

        @param n Anzahl der Elemente.
        @param k Anzahl der Anstiege (0 ≤ k < n).
        @return Euler-Zahl A(n,k).
        """
        # Randbedingungen
        if n == 0:
            return 1 if k == 0 else 0
        if k < 0 or k >= n:
            return 0

        # Rekurrenz für Euler-Zahlen
        return ((k + 1) * EnumerativeCombinatorics.euler_number_a(n - 1, k)
                + (n - k) * EnumerativeCombinatorics.euler_number_a(n - 1, k - 1))

    @staticmethod
    def partition_number(n: int) -> int:
        """
        @brief Berechnet die Partitionszahl p(n) – Anzahl der Partitionen von n.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        p(n) = Anzahl der Möglichkeiten, n als geordnete Summe positiver ganzer Zahlen
        zu schreiben (Reihenfolge irrelevant).

        Beispiel: p(4) = 5: 4 = 3+1 = 2+2 = 2+1+1 = 1+1+1+1.

        Methode: Dynamische Programmierung über Euler'sche Produktformel.
        p[n] = Σ_{k≠0} (-1)^{k+1} p[n - k(3k-1)/2]
        (Euler'sche Pentagonalzahl-Rekurrenz)

        @param n Nicht-negative ganze Zahl.
        @return Partitionszahl p(n).
        """
        if n < 0:
            return 0
        if n == 0:
            return 1

        # Tabelle für Partitionszahlen aufbauen (DP)
        p = [0] * (n + 1)
        p[0] = 1  # Basis: leere Partition

        # Euler'sche Pentagonalzahl-Rekurrenz
        for i in range(1, n + 1):
            total = 0
            k = 1
            while True:
                # Pentagonalzahlen: g_k = k(3k-1)/2, g_{-k} = k(3k+1)/2
                penta1 = k * (3 * k - 1) // 2
                penta2 = k * (3 * k + 1) // 2

                if penta1 > i:
                    break  # Alle relevanten Terme verarbeitet

                # Vorzeichen: (-1)^{k+1} → +1 für k=1,2; -1 für k=3,4; ...
                sign = 1 if k % 2 == 1 else -1
                total += sign * p[i - penta1]
                if penta2 <= i:
                    total += sign * p[i - penta2]

                k += 1

            p[i] = total

        return p[n]

    @staticmethod
    def motzkin_number(n: int) -> int:
        """
        @brief Berechnet die n-te Motzkin-Zahl M_n.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        M_n = Anzahl der Pfade von (0,0) nach (n,0) mit Schritten (+1,+1), (+1,0), (+1,-1),
        die nie unter die x-Achse gehen. Äquivalent: Anzahl nicht-kreuzender Partitionen
        mit Singleton-Teilen.

        Konvolutions-Rekurrenz: M_n = M_{n-1} + Σ_{k=0}^{n-2} M_k · M_{n-2-k}
        (Herleitung: jeder Motzkin-Pfad beginnt entweder flach (→ M_{n-1})
        oder steigt zuerst an, klettert k Schritte, und sinkt dann (→ M_k · M_{n-2-k}))

        M_0=1, M_1=1, M_2=2, M_3=4, M_4=9, M_5=21, ...

        @param n Index (n ≥ 0).
        @return n-te Motzkin-Zahl.
        """
        if n < 0:
            return 0
        if n <= 1:
            return 1

        # Iterative Berechnung via Konvolutionsrekurrenz (stabiler als rein rekursiv)
        M = [1, 1]  # M_0=1, M_1=1
        for i in range(2, n + 1):
            # M_i = M_{i-1} + Σ_{k=0}^{i-2} M_k * M_{i-2-k}
            conv_sum = sum(M[k] * M[i - 2 - k] for k in range(i - 1))
            M.append(M[i - 1] + conv_sum)

        return M[n]

    @staticmethod
    def narayana_number(n: int, k: int) -> int:
        """
        @brief Berechnet die Narayana-Zahl N(n,k).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        N(n,k) = (1/n) · C(n,k) · C(n,k-1)

        N(n,k) = Anzahl der nicht-kreuzenden Partitionen von {1,...,n} in k Blöcke.
        Auch: Anzahl der Dyck-Pfade der Länge 2n mit genau k Gipfeln (UD-Muster).

        Σ_{k=1}^{n} N(n,k) = C_n (n-te Catalan-Zahl).

        @param n Gesamtanzahl.
        @param k Anzahl der Blöcke/Gipfel (1 ≤ k ≤ n).
        @return Narayana-Zahl N(n,k).
        """
        if n <= 0 or k <= 0 or k > n:
            return 0
        if n == 1 and k == 1:
            return 1

        # N(n,k) = (1/n) * C(n,k) * C(n,k-1)
        return comb(n, k) * comb(n, k - 1) // n

    @staticmethod
    def ballot_problem(p: int, q: int) -> float:
        """
        @brief Berechnet die Wahrscheinlichkeit im Wahlproblem (Bertrand'sches Wahlproblem).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Kandidat A erhält p Stimmen, Kandidat B erhält q Stimmen (p > q).
        Gesucht: P(A liegt in jedem Zähltschritt vorn).

        Bertrand-Ballot-Theorem: P = (p - q) / (p + q).

        @param p Stimmen für Kandidat A (muss > q gelten).
        @param q Stimmen für Kandidat B.
        @return Wahrscheinlichkeit, dass A durchgehend in Führung liegt (0 wenn p ≤ q).
        """
        if p <= q:
            return 0.0  # A kann nicht immer vorne liegen, wenn er nicht gewinnt

        total = p + q
        # Bertrand-Ballot-Theorem: (p - q) / (p + q)
        return (p - q) / total

    @staticmethod
    def inclusion_exclusion(sets: list[set]) -> int:
        """
        @brief Berechnet |A_1 ∪ A_2 ∪ ... ∪ A_m| via Inklusion-Exklusion-Prinzip.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        |⋃ A_i| = Σ|A_i| - Σ|A_i ∩ A_j| + Σ|A_i ∩ A_j ∩ A_k| - ...

        Allgemein: Σ_{∅≠S⊆{1,...,m}} (-1)^{|S|+1} |⋂_{i∈S} A_i|

        @param sets Liste von Python-Mengen.
        @return Größe der Vereinigung.
        """
        m = len(sets)
        if m == 0:
            return 0

        total = 0
        # Über alle nicht-leeren Teilmengen iterieren
        for r in range(1, m + 1):
            for combo in itertools.combinations(range(m), r):
                # Schnittmenge aller Mengen in dieser Teilmenge
                intersection = sets[combo[0]].copy()
                for idx in combo[1:]:
                    intersection &= sets[idx]

                # Vorzeichen: (-1)^{r+1}
                sign = 1 if r % 2 == 1 else -1
                total += sign * len(intersection)

        return total


# ============================================================
# 3. GeneratingFunctions – Erzeugende Funktionen
# ============================================================

class GeneratingFunctions:
    """
    @brief Erzeugende Funktionen – gewöhnliche und exponentielle erzeugende Funktionen,
           Potenzreihenmultiplikation und Rekurrenzrelationen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Eine erzeugende Funktion ist eine formale Potenzreihe, deren Koeffizienten
    eine Folge kodieren.

    OGF: F(x) = Σ a_n x^n  (gewöhnliche erzeugende Funktion)
    EGF: F(x) = Σ a_n x^n/n!  (exponentielle erzeugende Funktion)
    """

    @staticmethod
    def ordinary_gf(sequence: list[float], n_terms: int) -> list[float]:
        """
        @brief Gibt die ersten n_terms Koeffizienten der gewöhnlichen erzeugenden Funktion.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        OGF(sequence) = a_0 + a_1·x + a_2·x² + ...
        Rückgabe: Liste der Koeffizienten [a_0, a_1, ..., a_{n_terms-1}].

        @param sequence Eingangsfolge.
        @param n_terms Anzahl der gewünschten Terme.
        @return Koeffizientenliste der Länge n_terms.
        """
        # Eingabe auffüllen oder kürzen auf n_terms Terme
        coeffs = list(sequence[:n_terms])
        # Fehlende Terme mit 0 auffüllen
        while len(coeffs) < n_terms:
            coeffs.append(0)
        return coeffs

    @staticmethod
    def exponential_gf(sequence: list[float], n_terms: int) -> list[float]:
        """
        @brief Gibt die ersten n_terms Koeffizienten der exponentiellen erzeugenden Funktion.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        EGF: F(x) = Σ a_n/n! · x^n
        Koeffizient von x^n in EGF = a_n / n!
        Rückgabe: [a_0/0!, a_1/1!, a_2/2!, ..., a_{n-1}/(n-1)!]

        @param sequence Eingangsfolge (a_n).
        @param n_terms Anzahl der gewünschten Terme.
        @return EGF-Koeffizientenliste der Länge n_terms.
        """
        coeffs = []
        for i in range(n_terms):
            if i < len(sequence):
                # EGF-Koeffizient: a_n / n!
                coeffs.append(sequence[i] / factorial(i))
            else:
                coeffs.append(0.0)
        return coeffs

    @staticmethod
    def fibonacci_gf() -> str:
        """
        @brief Gibt die geschlossene Form der erzeugenden Funktion der Fibonacci-Zahlen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        OGF der Fibonacci-Zahlen: F(x) = x / (1 - x - x²)

        Denn: F(x) - x·F(x) - x²·F(x) = x (nach F_0=0, F_1=1)
        ⟹ F(x) = x / (1 - x - x²)

        @return String-Darstellung der erzeugenden Funktion.
        """
        return "F(x) = x / (1 - x - x^2)"

    @staticmethod
    def catalan_gf() -> str:
        """
        @brief Gibt die geschlossene Form der erzeugenden Funktion der Catalan-Zahlen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        OGF der Catalan-Zahlen: C(x) = (1 - √(1 - 4x)) / (2x)

        Denn: C(x) = 1 + x·C(x)² (Dyck-Pfad-Rekurrenz)
        Quadratische Gleichung lösen: C(x) = (1 ± √(1-4x))/(2x)
        Reguläres Vorzeichen − ergibt die korrekte Lösung (C(0)=1).

        @return String-Darstellung der erzeugenden Funktion.
        """
        return "C(x) = (1 - sqrt(1 - 4x)) / (2x)"

    @staticmethod
    def power_series_multiply(f: list[float], g: list[float], n_terms: int) -> list[float]:
        """
        @brief Multipliziert zwei Potenzreihen (Cauchy-Faltung) und gibt n_terms Koeffizienten.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        (f · g)[n] = Σ_{k=0}^{n} f[k] · g[n-k]  (Cauchy-Faltung / Diskretes Konvolutionsprodukt)

        @param f Koeffizienten der ersten Potenzreihe.
        @param g Koeffizienten der zweiten Potenzreihe.
        @param n_terms Anzahl der gewünschten Ausgabeterme.
        @return Koeffizienten des Produkts f·g.
        """
        result = [0.0] * n_terms

        for n in range(n_terms):
            # Faltung: h[n] = Σ_{k=0}^{n} f[k] * g[n-k]
            for k in range(n + 1):
                f_k = f[k] if k < len(f) else 0.0
                g_nk = g[n - k] if (n - k) < len(g) else 0.0
                result[n] += f_k * g_nk

        return result

    @staticmethod
    def gf_from_recurrence(a0: float, a1: float, recurrence_fn, n_terms: int) -> list[float]:
        """
        @brief Berechnet eine Folge aus einer Rekurrenzrelation.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Startet mit a0, a1 und wendet die Rekurrenzfunktion an, um die Folge zu erzeugen.

        @param a0 Erstes Folgenglied a_0.
        @param a1 Zweites Folgenglied a_1.
        @param recurrence_fn Funktion f(seq, n) → a_n, erhält bisherige Folge und Index.
        @param n_terms Anzahl der zu berechnenden Terme.
        @return Folge der ersten n_terms Glieder.
        """
        if n_terms <= 0:
            return []
        if n_terms == 1:
            return [a0]

        seq = [a0, a1]

        # Rekurrenz iterativ anwenden
        for n in range(2, n_terms):
            next_val = recurrence_fn(seq, n)
            seq.append(next_val)

        return seq[:n_terms]


# ============================================================
# 4. RamseyTheory – Ramsey-Theorie
# ============================================================

class RamseyTheory:
    """
    @brief Ramsey-Theorie – Ramsey-Zahlen, Van-der-Waerden-Zahlen, Schubfachprinzip.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Ramsey-Satz: Für alle s,t ∈ ℕ existiert eine kleinste Zahl R(s,t), sodass
    jede 2-Färbung der Kanten von K_{R(s,t)} einen monochromatischen K_s oder K_t enthält.

    Bekannte Werte: R(3,3)=6, R(3,4)=9, R(3,5)=14, R(4,4)=18, ...
    """

    # Bekannte Ramsey-Zahlen R(s,t) für kleine s,t (s ≤ t)
    # Quelle: Radziszowski, "Small Ramsey Numbers", Electronic Journal of Combinatorics
    _RAMSEY_KNOWN: dict[tuple[int, int], int] = {
        (1, 1): 1, (1, 2): 2, (1, 3): 3, (1, 4): 4, (1, 5): 5, (1, 6): 6,
        (1, 7): 7, (1, 8): 8, (1, 9): 9, (1, 10): 10,
        (2, 2): 3, (2, 3): 3, (2, 4): 4, (2, 5): 5, (2, 6): 6,
        (2, 7): 7, (2, 8): 8, (2, 9): 9, (2, 10): 10,
        (3, 3): 6, (3, 4): 9, (3, 5): 14, (3, 6): 18, (3, 7): 23,
        (3, 8): 28, (3, 9): 36,
        (4, 4): 18, (4, 5): 25,
        (5, 5): 43,  # Untere Schranke, exakter Wert nicht bekannt
    }

    @classmethod
    def ramsey_number_R(cls, s: int, t: int) -> int | None:
        """
        @brief Gibt die Ramsey-Zahl R(s,t) zurück, falls bekannt.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        R(s,t) ist die kleinste Zahl n, sodass jede 2-Färbung der Kanten von K_n
        einen roten K_s oder blauen K_t enthält.

        Symmetrie: R(s,t) = R(t,s).
        Obere Schranke: R(s,t) ≤ R(s-1,t) + R(s,t-1) (Ramsey's Beweis via Induktion).

        @param s Größe des ersten monochromatischen Teilgraphen.
        @param t Größe des zweiten monochromatischen Teilgraphen.
        @return R(s,t) falls bekannt, None sonst.
        """
        # Symmetrie nutzen: s ≤ t
        if s > t:
            s, t = t, s

        # Triviale Fälle
        if s == 1:
            return t  # K_n enthält immer K_1

        # Bekannte Werte nachschlagen
        return cls._RAMSEY_KNOWN.get((s, t), None)

    @staticmethod
    def ramsey_graph_coloring(n: int, k: int) -> dict[tuple[int, int], int]:
        """
        @brief Erzeugt eine zufällige k-Färbung der Kanten des vollständigen Graphen K_n.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Jede Kante {u,v} des K_n erhält eine Farbe aus {0, 1, ..., k-1}.
        Diese deterministische Greedy-Färbung dient zur Demonstration.

        @param n Anzahl der Knoten.
        @param k Anzahl der Farben.
        @return Dictionary {(u,v): Farbe} für alle Kanten u < v.
        """
        coloring = {}
        color = 0

        # Alle Kanten des vollständigen Graphen K_n einfärben
        for u in range(n):
            for v in range(u + 1, n):
                # Farbe rotierend vergeben (Greedy-Zuweisung)
                coloring[(u, v)] = color % k
                color += 1

        return coloring

    @staticmethod
    def van_der_waerden_numbers() -> dict[tuple[int, int], int]:
        """
        @brief Gibt bekannte Van-der-Waerden-Zahlen W(k; r) zurück.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        W(k; r) = kleinste n, sodass jede r-Färbung von {1,...,n}
        eine monochromatische arithmetische Progression der Länge k enthält.

        Van-der-Waerden-Satz (1927): W(k; r) existiert für alle k, r.

        Notation: W(k; r) = W(k, k, ..., k) mit r Argumenten.
        Hier: W(k; 2) für 2 Farben.

        @return Dictionary {(k, r): W(k;r)} für bekannte Werte.
        """
        # Bekannte Van-der-Waerden-Zahlen (Quelle: OEIS, Landman & Robertson)
        return {
            (2, 2): 3,    # W(2;2)=3: {1,2,3} → immer AP der Länge 2
            (3, 2): 9,    # W(3;2)=9
            (4, 2): 35,   # W(4;2)=35
            (5, 2): 178,  # W(5;2)=178
            (6, 2): 1132, # W(6;2)=1132
            (3, 3): 27,   # W(3;3)=27
            (3, 4): 76,   # W(3;4)=76
        }

    @staticmethod
    def hales_jewett_theorem_demo(n: int, k: int) -> dict:
        """
        @brief Demonstriert das Hales-Jewett-Theorem für n-in-a-row auf k×k×...×k Hyperwürfel.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Hales-Jewett-Theorem (1963): Für beliebige k und t existiert d = HJ(k,t) so,
        dass jede k-Färbung des d-dimensionalen k-Gitters eine monochromatische
        kombinatorische Linie enthält.

        Dies impliziert Van-der-Waerden und verallgemeinert Tic-Tac-Toe.

        @param n Dimension.
        @param k Anzahl der Symbole/Farben.
        @return Demo-Information als Dictionary.
        """
        return {
            "theorem": "Hales-Jewett-Theorem (1963)",
            "dimension": n,
            "alphabet_size": k,
            "cube_size": k ** n,
            "description": (
                f"Jede {k}-Färbung des {n}-dimensionalen {k}-Gitters "
                f"({k}^{n} = {k**n} Punkte) enthält eine monochromatische "
                f"kombinatorische Linie."
            ),
            "implication": "Impliziert Van-der-Waerden-Theorem und Furstenberg-Katznelson.",
        }

    @staticmethod
    def pigeonhole_principle(n_pigeons: int, n_holes: int) -> dict:
        """
        @brief Schubfachprinzip: Wenn n_pigeons Tauben in n_holes Schubfächer verteilt werden.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Schubfachprinzip: Wenn n > k·m, dann muss mindestens ein Schubfach mehr als m Objekte enthalten.
        Einfache Form: n_pigeons > n_holes ⟹ mindestens ein Schubfach enthält ≥ 2 Tauben.

        Verallgemeinerung: Mindestens ein Schubfach enthält ≥ ⌈n_pigeons/n_holes⌉ Tauben.

        @param n_pigeons Anzahl der Tauben (Objekte).
        @param n_holes Anzahl der Schubfächer (Kategorien).
        @return Dictionary mit Analyse und Garantien.
        """
        import math as _math
        if n_holes <= 0:
            raise ValueError("Anzahl der Schubfächer muss positiv sein.")

        # Mindestzahl Tauben in einem Schubfach (Deckenftunktion)
        min_max_pigeons = _math.ceil(n_pigeons / n_holes)

        return {
            "n_pigeons": n_pigeons,
            "n_holes": n_holes,
            "guaranteed_min_in_one_hole": min_max_pigeons,
            "at_least_two_in_one_hole": n_pigeons > n_holes,
            "principle": (
                f"Bei {n_pigeons} Tauben und {n_holes} Schubfächern "
                f"enthält mindestens ein Schubfach ≥ {min_max_pigeons} Tauben."
            ),
        }


# ============================================================
# 5. GraphCombinatorics – Graphen-Kombinatorik
# ============================================================

class GraphCombinatorics:
    """
    @brief Graphen-Kombinatorik – chromatische Polynome, Matchings, Spannbäume,
           Euler- und Hamilton-Kreise.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Adjazenzmatrix-Darstellung: adj[i][j] = 1 wenn Kante {i,j} vorhanden, sonst 0.
    Ungerichtete Graphen: adj[i][j] = adj[j][i].
    """

    @staticmethod
    def chromatic_polynomial(adj_matrix: list[list[int]]) -> list[int]:
        """
        @brief Berechnet das chromatische Polynom P(G, k) via Deletion-Contraction.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        P(G, k) = Anzahl der k-Färbungen von G (ordnungsgemäße Knotenfärbungen).
        P(G, k) ist ein Polynom in k vom Grad n = |V(G)|.

        Deletion-Contraction: Für jede Kante e = {u,v}:
            P(G, k) = P(G - e, k) - P(G / e, k)
        Dabei ist G-e der Graph ohne Kante e, G/e der kontrahierte Graph.

        Basisfall: Leerer Graph K_n^0 hat P = k^n (jeder Knoten frei wählbar).

        Rückgabe: Koeffizientenliste [a_0, a_1, ..., a_n] mit P(k) = Σ a_i k^i.

        @param adj_matrix Quadratische 0/1-Adjazenzmatrix.
        @return Koeffizienten des chromatischen Polynoms.
        """
        n = len(adj_matrix)
        if n == 0:
            # Leerer Graph: P = 1 (das "leere Produkt")
            return [1]

        # Kanten aus der Adjazenzmatrix extrahieren (nur i < j, ungerichtet)
        edges = []
        for i in range(n):
            for j in range(i + 1, n):
                if adj_matrix[i][j]:
                    edges.append((i, j))

        # Rekursive Berechnung via Deletion-Contraction
        def poly_add(p1: list[int], p2: list[int]) -> list[int]:
            """Addiert zwei Polynome (Koeffizientenlisten)."""
            size = max(len(p1), len(p2))
            result = [0] * size
            for idx, c in enumerate(p1):
                result[idx] += c
            for idx, c in enumerate(p2):
                result[idx] += c
            return result

        def poly_sub(p1: list[int], p2: list[int]) -> list[int]:
            """Subtrahiert zwei Polynome."""
            size = max(len(p1), len(p2))
            result = [0] * size
            for idx, c in enumerate(p1):
                result[idx] += c
            for idx, c in enumerate(p2):
                result[idx] -= c
            return result

        def chromatic_rec(n_nodes: int, edge_set: list[tuple[int, int]]) -> list[int]:
            """Rekursive Deletion-Contraction Berechnung."""
            # Kein Kantengraph: P = k^n
            if not edge_set:
                # k^n als Polynom: [0, 0, ..., 0, 1] mit 1 an Position n_nodes
                poly = [0] * (n_nodes + 1)
                poly[n_nodes] = 1
                return poly

            # Erste Kante wählen
            u, v = edge_set[0]
            remaining = edge_set[1:]

            # Deletion: G - e (Kante entfernen)
            p_deletion = chromatic_rec(n_nodes, remaining)

            # Contraction: G / e (Knoten u und v verschmelzen)
            # Knoten umnummerieren: v wird zu u, alle v-Kanten werden zu u-Kanten
            contracted_edges = []
            seen = set()
            for a, b in remaining:
                # v durch u ersetzen
                a2 = u if a == v else a
                b2 = u if b == v else b
                if a2 != b2:  # Schleifen entfernen
                    edge_key = (min(a2, b2), max(a2, b2))
                    if edge_key not in seen:
                        seen.add(edge_key)
                        contracted_edges.append(edge_key)

            p_contraction = chromatic_rec(n_nodes - 1, contracted_edges)

            # Deletion-Contraction: P(G) = P(G-e) - P(G/e)
            return poly_sub(p_deletion, p_contraction)

        return chromatic_rec(n, edges)

    @staticmethod
    def tutte_polynomial_demo(adj_matrix: list[list[int]]) -> dict:
        """
        @brief Demo des Tutte-Polynoms T(G; x, y) für einfache Graphen.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Das Tutte-Polynom ist eine 2-variable Verallgemeinerung des chromatischen Polynoms.
        Spezielle Auswertungen:
          - T(G; 1, 1) = Anzahl der Spannbäume
          - T(G; 2, 2) = 2^|E| (Anzahl aller Teilgraphen)
          - T(G; 0, 2) = Anzahl der Zyklen (orientiert)
          - P(G, k) = (-1)^{n-c} k^c T(G; 1-k, 0)

        Diese Demo berechnet T(G; 1,1) = Spannbaumanzahl.

        @param adj_matrix Quadratische 0/1-Adjazenzmatrix.
        @return Demo-Dictionary mit Spezialwerten.
        """
        n = len(adj_matrix)
        spanning_trees = GraphCombinatorics.spanning_tree_count(adj_matrix)
        edges = sum(1 for i in range(n) for j in range(i + 1, n) if adj_matrix[i][j])

        return {
            "n_vertices": n,
            "n_edges": edges,
            "T_1_1": spanning_trees,  # T(G; 1,1) = Spannbaumanzahl
            "T_2_2": 2 ** edges,      # T(G; 2,2) = 2^|E|
            "description": "Tutte-Polynom T(G; x, y): T(1,1) = Spannbäume, T(2,2) = 2^|E|",
        }

    @staticmethod
    def matching_number(adj_matrix: list[list[int]]) -> int:
        """
        @brief Berechnet die Matchingzahl ν(G) – Größe des maximalen Matchings.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein Matching M ⊆ E ist eine Menge knotendisjunkter Kanten.
        ν(G) = max{|M| : M ist Matching in G}.

        Methode: Augmentierende-Pfade-Algorithmus (vereinfachte Version via Greedy + Augmentation).
        Für bipartite Graphen: exakt via Hopcroft-Karp (hier Greedy als Näherung).

        @param adj_matrix Quadratische 0/1-Adjazenzmatrix.
        @return Größe des maximalen Matchings.
        """
        n = len(adj_matrix)
        if n == 0:
            return 0

        matched = [-1] * n  # matched[v] = Partner von v, oder -1

        def try_augment(v: int, visited: list[bool]) -> bool:
            """Versucht, augmentierenden Pfad von v zu finden (DFS)."""
            for u in range(n):
                if adj_matrix[v][u] and not visited[u]:
                    visited[u] = True
                    # Freier Knoten oder Augmentation möglich
                    if matched[u] == -1 or try_augment(matched[u], visited):
                        matched[u] = v
                        return True
            return False

        # Maximales Matching via augmentierende Pfade
        match_count = 0
        for v in range(n):
            visited = [False] * n
            if try_augment(v, visited):
                match_count += 1

        # Jede Kante wurde doppelt gezählt (u→v und v→u)
        return match_count // 2

    @staticmethod
    def perfect_matching_count(adj_matrix: list[list[int]]) -> int:
        """
        @brief Berechnet die Anzahl perfekter Matchings via Hafnian (allgemeine Graphen).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein perfektes Matching M umfasst alle n Knoten (n muss gerade sein) als
        paarweise disjunkte Kanten.

        Für allgemeine ungerichtete Graphen wird der Hafnian verwendet:
            haf(A) = Σ_{M perfektes Matching} Π_{(i,j)∈M} A_{i,j}

        Methode: Rekursive Aufzählung aller perfekten Matchings (für kleine n).
        Komplexität: O(n!! · n) – in der Praxis exponentiell.

        Hinweis: Die Permanente gilt für bipartite Graphen. Für allgemeine Graphen
        ist der Hafnian das korrekte Maß.

        @param adj_matrix Quadratische 0/1-Adjazenzmatrix (n×n, n gerade für perfektes Matching).
        @return Anzahl der perfekten Matchings.
        """
        n = len(adj_matrix)
        if n == 0:
            return 1
        if n % 2 != 0:
            return 0  # Perfektes Matching nur bei gerader Knotenanzahl möglich

        def count_matchings(available: list[int]) -> int:
            """Rekursive Aufzählung perfekter Matchings via Backtracking."""
            if not available:
                return 1  # Leere Knotenmenge: ein (leeres) Matching gefunden

            # Kleinsten verfügbaren Knoten als Ankerpunkt wählen
            first = available[0]
            rest = available[1:]
            total = 0

            # Alle möglichen Partner für 'first' durchprobieren
            for i, partner in enumerate(rest):
                if adj_matrix[first][partner]:
                    # Kante {first, partner} existiert → in Matching aufnehmen
                    new_available = [x for x in rest if x != partner]
                    total += count_matchings(new_available)

            return total

        # Alle Knotenindizes als Ausgangsmenge
        nodes = list(range(n))
        return count_matchings(nodes)

    @staticmethod
    def spanning_tree_count(adj_matrix: list[list[int]]) -> int:
        """
        @brief Berechnet die Anzahl der Spannbäume via Kirchhoff-Matrix-Satz.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Kirchhoff-Matrix-Satz (Matrix-Baum-Satz):
        τ(G) = det(L̃)  wobei L̃ die reduzierte Laplace-Matrix ist
        (L ohne eine beliebige Zeile und Spalte).

        Laplace-Matrix: L = D - A
          D: Gradmatrix (d_{ii} = Grad von Knoten i)
          A: Adjazenzmatrix

        @param adj_matrix Quadratische 0/1-Adjazenzmatrix.
        @return Anzahl der Spannbäume (0 wenn nicht zusammenhängend).
        """
        n = len(adj_matrix)
        if n == 0:
            return 0
        if n == 1:
            return 1  # Einzelknoten hat genau 1 Spannbaum (er selbst)

        # Laplace-Matrix berechnen: L = D - A
        laplace = [[0] * n for _ in range(n)]
        for i in range(n):
            degree = sum(adj_matrix[i])
            laplace[i][i] = degree  # Diagonale: Grad des Knotens
            for j in range(n):
                if adj_matrix[i][j]:
                    laplace[i][j] -= 1  # Off-Diagonal: -1 wenn Kante vorhanden

        # Reduzierte Laplace-Matrix: letzte Zeile und Spalte entfernen
        reduced = [[laplace[i][j] for j in range(n - 1)] for i in range(n - 1)]

        # Determinante der reduzierten Matrix berechnen (Gauss-Elimination)
        def determinant_int(matrix: list[list[float]]) -> float:
            """Berechnet Determinante via Gauss-Elimination."""
            m = len(matrix)
            mat = [row[:] for row in matrix]  # Kopie
            det = 1.0

            for col in range(m):
                # Pivotzeile suchen (maximaler Betrag)
                pivot_row = max(range(col, m), key=lambda r: abs(mat[r][col]))
                if pivot_row != col:
                    mat[col], mat[pivot_row] = mat[pivot_row], mat[col]
                    det *= -1  # Zeilentausch ändert Vorzeichen

                pivot = mat[col][col]
                if abs(pivot) < 1e-10:
                    return 0.0  # Singuläre Matrix → det = 0

                det *= pivot

                # Elimination nach unten
                for row in range(col + 1, m):
                    factor = mat[row][col] / pivot
                    for k in range(col, m):
                        mat[row][k] -= factor * mat[col][k]

            return det

        # Determinante berechnen und auf ganze Zahl runden
        det = determinant_int([[float(reduced[i][j]) for j in range(n - 1)]
                               for i in range(n - 1)])
        return round(abs(det))

    @staticmethod
    def eulerian_circuit_exists(adj_matrix: list[list[int]]) -> bool:
        """
        @brief Prüft, ob der Graph einen Euler-Kreis besitzt.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Euler-Kreis-Satz (Euler 1736): Ein zusammenhängender Graph hat genau dann
        einen Euler-Kreis, wenn alle Knoten geraden Grad haben.

        Methode: Überprüfe (1) Zusammenhang und (2) alle Grade gerade.

        @param adj_matrix Quadratische 0/1-Adjazenzmatrix.
        @return True wenn Euler-Kreis existiert, False sonst.
        """
        n = len(adj_matrix)
        if n == 0:
            return True  # Leerer Graph: trivial

        # Bedingung 1: Alle Knoten müssen geraden Grad haben
        for i in range(n):
            degree = sum(adj_matrix[i])
            if degree % 2 != 0:
                return False  # Ungerader Grad → kein Euler-Kreis

        # Bedingung 2: Graph muss zusammenhängend sein (nur Knoten mit Grad > 0 betrachten)
        start = next((i for i in range(n) if sum(adj_matrix[i]) > 0), None)
        if start is None:
            return True  # Kein Knoten mit Kanten → trivial

        # BFS/DFS Zusammenhangs-Check
        visited = [False] * n
        stack = [start]
        visited[start] = True
        while stack:
            v = stack.pop()
            for u in range(n):
                if adj_matrix[v][u] and not visited[u]:
                    visited[u] = True
                    stack.append(u)

        # Alle Knoten mit Grad > 0 müssen besucht worden sein
        for i in range(n):
            if sum(adj_matrix[i]) > 0 and not visited[i]:
                return False  # Nicht zusammenhängend

        return True

    @staticmethod
    def hamiltonian_path_check(adj_matrix: list[list[int]]) -> bool:
        """
        @brief Prüft via Backtracking, ob der Graph einen Hamiltonpfad enthält.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein Hamiltonpfad besucht jeden Knoten genau einmal.
        Das Hamiltonpfad-Problem ist NP-vollständig.

        Methode: Backtracking mit Pruning (exponentielle Worst-Case-Komplexität O(n!)).

        @param adj_matrix Quadratische 0/1-Adjazenzmatrix.
        @return True wenn Hamiltonpfad existiert, False sonst.
        """
        n = len(adj_matrix)
        if n == 0:
            return True
        if n == 1:
            return True  # Einzelknoten ist selbst ein Hamiltonpfad

        visited = [False] * n

        def backtrack(current: int, count: int) -> bool:
            """Rekursive Backtracking-Suche."""
            if count == n:
                return True  # Alle Knoten besucht → Hamiltonpfad gefunden

            for next_v in range(n):
                if adj_matrix[current][next_v] and not visited[next_v]:
                    visited[next_v] = True
                    if backtrack(next_v, count + 1):
                        return True
                    visited[next_v] = False  # Backtrack

            return False  # Kein Hamiltonpfad von current aus

        # Von jedem Startknoten versuchen
        for start in range(n):
            visited[start] = True
            if backtrack(start, 1):
                return True
            visited[start] = False

        return False


# ============================================================
# 6. CombinatoricsOnWords – Kombinatorik auf Wörtern
# ============================================================

class CombinatoricsOnWords:
    """
    @brief Kombinatorik auf Wörtern – Lyndon-Wörter, Halsketten, LCS, Edit-Distanz, De-Bruijn.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Ein Wort w über Alphabet Σ ist ein Element Σ*.
    Periodizität: w = u^k bedeutet w ist k-fache Wiederholung von u.
    Primitiv: w ist primitiv, wenn kein echtes Teilwort u existiert mit w = u^k.
    Lyndon-Wort: primitives Wort, das lexikographisch kleiner ist als alle seine Rotationen.
    """

    @staticmethod
    def lyndon_word_check(word: str) -> bool:
        """
        @brief Prüft, ob ein Wort ein Lyndon-Wort ist.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein Lyndon-Wort ist ein nicht-leeres Wort w, das:
        1. Primitiv ist (keine echte Periode).
        2. Lexikographisch kleiner ist als alle seine echten, nicht-leeren Rotationen.

        Methode: Vergleiche w mit allen Rotationen rot_i(w) = w[i:] + w[:i].

        @param word Eingabewort als String.
        @return True wenn Lyndon-Wort, False sonst.
        """
        if not word:
            return False  # Leeres Wort ist kein Lyndon-Wort per Definition

        n = len(word)

        # Alle echten Rotationen prüfen
        for i in range(1, n):
            rotation = word[i:] + word[:i]
            if rotation <= word:
                # Rotation ist kleiner oder gleich → kein Lyndon-Wort
                return False

        return True  # Alle Rotationen sind lexikographisch größer

    @staticmethod
    def necklaces(n: int, k: int) -> int:
        """
        @brief Berechnet die Anzahl der Halsketten mit n Perlen aus k Farben (Burnside-Lemma).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Eine Halskette ist eine Äquivalenzklasse von Wörtern der Länge n unter Rotation.
        Burnside-Lemma: |X/G| = (1/|G|) Σ_{g∈G} |Fix(g)|

        Für zyklische Gruppe Z_n: |Fix(Rotation um d)| = k^{ggT(n,d)}.
        Anzahl Halsketten = (1/n) Σ_{d=0}^{n-1} k^{ggT(n,d)}.

        @param n Anzahl der Perlen (Länge des Wortes).
        @param k Anzahl der Farben/Symbole.
        @return Anzahl der verschiedenen Halsketten.
        """
        if n == 0:
            return 1  # Leere Halskette
        if k == 0:
            return 0  # Keine Farben → keine Halsketten

        # Burnside-Lemma: (1/n) Σ_{d=0}^{n-1} k^{ggT(n,d)}
        total = sum(k ** gcd(n, d) for d in range(n))
        return total // n

    @staticmethod
    def primitive_word(word: str) -> bool:
        """
        @brief Prüft, ob ein Wort primitiv ist (keine echte Periode).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Ein Wort w ist primitiv, wenn kein echtes Teilwort u ≠ w und k ≥ 2 existiert
        mit w = u^k (w also keine Wiederholung eines kürzeren Wortes ist).

        Methode: Prüfe alle echten Teiler der Länge n auf Periodizität.

        @param word Eingabewort.
        @return True wenn primitiv, False wenn w = u^k für k ≥ 2.
        """
        n = len(word)
        if n <= 1:
            return True  # Einzelzeichen sind immer primitiv

        # Alle echten Teiler von n prüfen
        for period_len in range(1, n):
            if n % period_len == 0:
                # Prüfe ob word = word[:period_len] ^ (n/period_len)
                period = word[:period_len]
                k = n // period_len
                if period * k == word:
                    if period_len < n:
                        return False  # Nicht primitiv

        return True  # Kein echtes Periode-Muster gefunden

    @staticmethod
    def longest_common_subsequence(s1: str, s2: str) -> int:
        """
        @brief Berechnet die Länge der längsten gemeinsamen Teilfolge (LCS).
        @author Kurt Ingwer
        @lastModified 2026-03-10

        LCS(s1, s2) = längste Folge, die in s1 und s2 als Teilfolge (nicht zwingend zusammenhängend)
        vorkommt.

        Methode: Dynamische Programmierung in O(m·n) Zeit und Speicher.
        Rekurrenz:
          dp[i][j] = dp[i-1][j-1] + 1       falls s1[i] = s2[j]
          dp[i][j] = max(dp[i-1][j], dp[i][j-1])  sonst

        @param s1 Erste Zeichenkette.
        @param s2 Zweite Zeichenkette.
        @return Länge der LCS.
        """
        m, n = len(s1), len(s2)

        # DP-Tabelle initialisieren (m+1) × (n+1)
        dp = [[0] * (n + 1) for _ in range(m + 1)]

        # DP-Tabelle füllen
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if s1[i - 1] == s2[j - 1]:
                    # Zeichen stimmt überein: diagonal verlängern
                    dp[i][j] = dp[i - 1][j - 1] + 1
                else:
                    # Kein Match: Maximum der Nachbarn
                    dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

        return dp[m][n]

    @staticmethod
    def edit_distance(s1: str, s2: str) -> int:
        """
        @brief Berechnet die Levenshtein-Distanz (Edit-Distanz) zwischen zwei Strings.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Die Edit-Distanz ist die minimale Anzahl von Einfügungen, Löschungen und
        Ersetzungen, um s1 in s2 umzuwandeln.

        Methode: Wagner-Fischer-Algorithmus (DP) in O(m·n).
        Rekurrenz:
          dp[i][j] = 0                            falls i=0 oder j=0 Basis
          dp[i][j] = dp[i-1][j-1]                falls s1[i]=s2[j]
          dp[i][j] = 1 + min(dp[i-1][j],         (Löschung)
                              dp[i][j-1],         (Einfügung)
                              dp[i-1][j-1])       (Ersetzung)

        @param s1 Erste Zeichenkette.
        @param s2 Zweite Zeichenkette.
        @return Levenshtein-Distanz.
        """
        m, n = len(s1), len(s2)

        # Spezialfälle
        if m == 0:
            return n  # s1 leer → n Einfügungen nötig
        if n == 0:
            return m  # s2 leer → m Löschungen nötig

        # DP-Tabelle: dp[i][j] = Edit-Distanz zwischen s1[:i] und s2[:j]
        dp = [[0] * (n + 1) for _ in range(m + 1)]

        # Basisfall: Transformationen leerer Strings
        for i in range(m + 1):
            dp[i][0] = i  # s2 leer: i Löschungen
        for j in range(n + 1):
            dp[0][j] = j  # s1 leer: j Einfügungen

        # Tabelle füllen
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if s1[i - 1] == s2[j - 1]:
                    # Zeichen identisch: keine Operation nötig
                    dp[i][j] = dp[i - 1][j - 1]
                else:
                    # Minimum aus Löschung, Einfügung, Ersetzung
                    dp[i][j] = 1 + min(
                        dp[i - 1][j],    # Löschung aus s1
                        dp[i][j - 1],    # Einfügung in s1
                        dp[i - 1][j - 1] # Ersetzung
                    )

        return dp[m][n]

    @staticmethod
    def de_bruijn_sequence(n: int, k: int) -> str:
        """
        @brief Generiert eine De-Bruijn-Folge B(k, n) – enthält jedes k-mere der Länge n genau einmal.
        @author Kurt Ingwer
        @lastModified 2026-03-10

        Eine De-Bruijn-Folge der Ordnung n über einem k-elementigen Alphabet ist eine
        zyklische Folge, in der jede mögliche Teilfolge der Länge n genau einmal vorkommt.
        Gesamtlänge: k^n Zeichen.

        Methode: Martin'scher Algorithmus (via Lynch-Rote-Rückruf-Algorithmus).
        Basiert auf dem eulerschen Kreis im De-Bruijn-Graphen.

        Alphabet: '0', '1', ..., str(k-1) für k ≤ 10; für k > 10 als Ziffernliste.

        @param n Länge der k-mere (Ordnung).
        @param k Alphabetgröße.
        @return De-Bruijn-Folge als String.
        """
        # Martin'scher De-Bruijn-Algorithmus
        alphabet = [str(i) for i in range(k)]

        # Resultat-Zeichenliste
        sequence = []
        a = [0] * (k * n)

        def db(t: int, p: int) -> None:
            """Rekursiver Lyndon-Wort-Aufbau (Martin-Algorithmus)."""
            if t > n:
                if n % p == 0:
                    # Gültiger Abschnitt: a[1..p] zur Folge hinzufügen
                    for i in range(1, p + 1):
                        sequence.append(alphabet[a[i]])
            else:
                a[t] = a[t - p]
                db(t + 1, p)
                for j in range(a[t - p] + 1, k):
                    a[t] = j
                    db(t + 1, t)

        db(1, 1)
        return ''.join(sequence)


# ============================================================
# 7. Standalone-Funktionen – Zahlenfolgen
# ============================================================

def fibonacci(n: int) -> int:
    """
    @brief Berechnet die n-te Fibonacci-Zahl via Matrix-Exponentiation in O(log n).
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Fibonacci-Folge: F_0=0, F_1=1, F_n=F_{n-1}+F_{n-2}.

    Matrix-Exponentiation: [[F_{n+1}], [F_n]] = [[1,1],[1,0]]^n · [[1],[0]]
    Dadurch O(log n) statt O(n) Multiplikationen.

    @param n Index der Fibonacci-Zahl (n ≥ 0).
    @return n-te Fibonacci-Zahl F_n.
    @raises ValueError für n < 0.
    """
    if n < 0:
        raise ValueError(f"Fibonacci ist nur für n ≥ 0 definiert, nicht für n={n}.")

    # Basisfall
    if n == 0:
        return 0
    if n == 1:
        return 1

    def mat_mul(A: list[list[int]], B: list[list[int]]) -> list[list[int]]:
        """2×2 Matrizen-Multiplikation."""
        return [
            [A[0][0] * B[0][0] + A[0][1] * B[1][0],
             A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0],
             A[1][0] * B[0][1] + A[1][1] * B[1][1]],
        ]

    def mat_pow(M: list[list[int]], exp: int) -> list[list[int]]:
        """Schnelle Matrizenpotenz via Binäre Exponentiation."""
        # Einheitsmatrix als Startwert
        result = [[1, 0], [0, 1]]
        base = M

        while exp > 0:
            if exp % 2 == 1:
                result = mat_mul(result, base)
            base = mat_mul(base, base)
            exp //= 2

        return result

    # Fibonacci-Matrix [[1,1],[1,0]]^n berechnen
    fib_matrix = [[1, 1], [1, 0]]
    result_matrix = mat_pow(fib_matrix, n)

    # F_n = Matrixelement [0][1] (oder [1][0])
    return result_matrix[0][1]


def lucas_numbers(n: int) -> list[int]:
    """
    @brief Berechnet die ersten n Lucas-Zahlen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Lucas-Folge: L_0=2, L_1=1, L_n=L_{n-1}+L_{n-2}.
    Verwandt mit Fibonacci: L_n = F_{n-1} + F_{n+1}.

    @param n Anzahl der gewünschten Lucas-Zahlen.
    @return Liste [L_0, L_1, ..., L_{n-1}].
    """
    if n <= 0:
        return []
    if n == 1:
        return [2]

    result = [2, 1]  # L_0=2, L_1=1

    # Restliche Terme via Rekurrenz
    for _ in range(2, n):
        result.append(result[-1] + result[-2])

    return result[:n]


def tribonacci(n: int) -> list[int]:
    """
    @brief Berechnet die ersten n Tribonacci-Zahlen.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Tribonacci-Folge: T_0=0, T_1=0, T_2=1, T_n=T_{n-1}+T_{n-2}+T_{n-3}.
    Verallgemeinert Fibonacci auf drei Vorgänger.

    @param n Anzahl der gewünschten Tribonacci-Zahlen.
    @return Liste [T_0, ..., T_{n-1}].
    """
    if n <= 0:
        return []
    if n == 1:
        return [0]
    if n == 2:
        return [0, 0]
    if n == 3:
        return [0, 0, 1]

    result = [0, 0, 1]

    # Rekurrenz: T_n = T_{n-1} + T_{n-2} + T_{n-3}
    for _ in range(3, n):
        result.append(result[-1] + result[-2] + result[-3])

    return result[:n]


def sylvester_sequence(n: int) -> list[int]:
    """
    @brief Berechnet die ersten n Glieder der Sylvester-Folge.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Sylvester-Folge: a_0=2, a_n=a_{n-1}·(a_{n-1}-1)+1 = a_0·a_1·...·a_{n-1}+1.
    Wächst doppelt-exponentiell.

    Verwendung in ägyptischen Bruchdarstellungen:
    1 = 1/a_0 + 1/(a_0·a_1) + ... (Sylvester-Entwicklung)

    a_0=2, a_1=3, a_2=7, a_3=43, a_4=1807, ...

    @param n Anzahl der gewünschten Glieder.
    @return Liste [a_0, ..., a_{n-1}].
    """
    if n <= 0:
        return []

    result = [2]  # a_0 = 2

    # Rekurrenz: a_n = a_{n-1} * (a_{n-1} - 1) + 1
    for _ in range(1, n):
        last = result[-1]
        result.append(last * (last - 1) + 1)

    return result[:n]


def happy_numbers(n: int) -> list[int]:
    """
    @brief Findet alle Glückszahlen bis n.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Eine positive ganze Zahl ist glücklich (happy), wenn die iterierte Anwendung
    der Quersumme der Quadratziffern schließlich zu 1 führt:
    z.B. 13 → 1²+3²=10 → 1²+0²=1 ✓

    Unglückliche Zahlen enden in einem Zyklus, der 4 enthält.

    @param n Obere Grenze (inklusive).
    @return Sortierte Liste aller Glückszahlen von 1 bis n.
    """
    def is_happy(num: int) -> bool:
        """Prüft ob eine Zahl glücklich ist (Floyd's Cycle Detection)."""
        seen = set()
        while num != 1 and num not in seen:
            seen.add(num)
            # Summe der Quadratziffern
            num = sum(int(d) ** 2 for d in str(num))
        return num == 1

    # Alle Glückszahlen von 1 bis n finden
    return [i for i in range(1, n + 1) if is_happy(i)]


def collatz_tree(n_max: int) -> dict[int, list[int]]:
    """
    @brief Erstellt den Collatz-Baum als Dictionary für alle Startpunkte bis n_max.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Der Collatz-Baum zeigt, welche Zahlen direkt zu einem Knoten führen
    (Vorgängerrelation in der Collatz-Folge).

    Collatz-Abbildung: n → n/2 (n gerade), n → 3n+1 (n ungerade).
    Umgekehrt: n ← 2n (immer), n ← (n-1)/3 (wenn n ≡ 1 mod 3 und (n-1)/3 ungerade).

    @param n_max Maximaler Startwert.
    @return Dictionary {n: [Vorgänger von n im Collatz-Graphen]}.
    """
    # Collatz-Baum als Vorgänger-Dictionary aufbauen
    tree: dict[int, list[int]] = {i: [] for i in range(1, n_max + 1)}

    for m in range(1, n_max + 1):
        # Nachfolger von m berechnen
        if m % 2 == 0:
            successor = m // 2
        else:
            successor = 3 * m + 1

        # m ist Vorgänger von successor
        if successor in tree:
            tree[successor].append(m)

    return tree
