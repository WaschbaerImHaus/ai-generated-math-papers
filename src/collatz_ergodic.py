"""
@file collatz_ergodic.py
@brief Collatz-Funktion, 2-adische Analyse und ergodische Maße.
@description
    Implementiert die Collatz-Vermutung aus verschiedenen mathematischen
    Blickwinkeln:

    1. **Klassisch**: Die Collatz-Funktion T: ℕ → ℕ mit
       T(n) = n/2 falls n gerade, 3n+1 falls n ungerade.
       Vermutung: Für alle n ≥ 1 erreicht die Folge T^k(n) schließlich 1.

    2. **2-adisch (nach Lagarias)**: Collatz als Funktion auf den 2-adischen
       Zahlen ℤ₂. Der Lyapunov-Exponent
           λ = log(3)/2 − log(2)/2 ≈ −0.0959 < 0
       erklärt die typische Konvergenz der Folge.

    3. **Ergodisch (nach Tao, 2019)**: Terenz Tao bewies 2019, dass fast alle
       Startwerte n eine Trajektorie haben, die beliebig oft unter jede
       divergierende Funktion f(n) fällt. Insbesondere gilt für jedes ε > 0:
           Dichte({n : min_k T^k(n) ≤ n^ε}) = 1.

    Mathematische Hintergründe:
    - Lagarias, J.C. (1985): "The 3x+1 problem and its generalizations"
    - Tao, T. (2019): "Almost all orbits of the Collatz map attain almost
      bounded values", arXiv:1909.03562

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import random
import numpy as np
import matplotlib.pyplot as plt
from typing import Callable


# ===========================================================================
# Klasse: CollatzFunction
# ===========================================================================

class CollatzFunction:
    """
    Grundlegende Collatz-Operationen.

    Die Collatz-Funktion ist definiert als:
        T(n) = n/2       falls n gerade
        T(n) = 3n + 1    falls n ungerade

    Die Collatz-Vermutung besagt: Für alle n ≥ 1 existiert ein k mit T^k(n) = 1.

    @author Michael Fuhrmann
    @since 2026-03-11
    """

    def __init__(self) -> None:
        """
        Initialisiert die Collatz-Funktion.
        Kein interner Zustand nötig – alle Methoden sind zustandslos.
        """
        pass

    def step(self, n: int) -> int:
        """
        Führt einen einzelnen Collatz-Schritt aus.

        Formel:
            T(n) = n/2       falls 2 | n
            T(n) = 3n + 1    sonst

        @param n: Positive ganze Zahl (Eingabe)
        @return: T(n)
        @raises ValueError: Falls n ≤ 0
        """
        if n <= 0:
            raise ValueError(f"Collatz ist nur für positive Zahlen definiert, erhielt: {n}")
        # Gerade: halbieren
        if n % 2 == 0:
            return n // 2
        # Ungerade: 3n+1
        return 3 * n + 1

    def sequence(self, n: int, max_steps: int = 10000) -> list:
        """
        Berechnet die vollständige Collatz-Folge ab n bis 1 (oder max_steps).

        Die Folge lautet: n, T(n), T²(n), T³(n), ..., 1

        @param n:         Startwert (positive ganze Zahl)
        @param max_steps: Maximale Anzahl von Schritten (Sicherheitsabbruch)
        @return:          Liste der Folgenglieder, inklusive Startwert und 1
        @raises ValueError: Falls n ≤ 0
        """
        if n <= 0:
            raise ValueError(f"Startwert muss positiv sein, erhielt: {n}")

        folge = [n]
        aktuell = n
        schritte = 0

        # Iteration bis 1 oder Schrittlimit
        while aktuell != 1 and schritte < max_steps:
            aktuell = self.step(aktuell)
            folge.append(aktuell)
            schritte += 1

        return folge

    def stopping_time(self, n: int) -> int:
        """
        Berechnet die Stoppzeit σ(n): Anzahl Schritte bis T^k(n) < n.

        Die Stoppzeit ist die kleinste Zahl k mit T^k(n) < n.
        Für n = 1: σ(1) = 0 (bereits minimal).

        @param n: Startwert (positive ganze Zahl)
        @return:  Stoppzeit σ(n); -1 falls kein k in 10000 Schritten gefunden
        """
        if n <= 1:
            return 0

        aktuell = n
        for k in range(1, 10001):
            aktuell = self.step(aktuell)
            if aktuell < n:
                return k

        return -1  # Nicht konvergiert (praktisch nie)

    def total_stopping_time(self, n: int) -> int:
        """
        Berechnet die totale Stoppzeit τ(n): Anzahl Schritte bis T^k(n) = 1.

        @param n: Startwert (positive ganze Zahl)
        @return:  Totale Stoppzeit τ(n)
        @raises ValueError: Falls n ≤ 0
        """
        if n <= 0:
            raise ValueError(f"Startwert muss positiv sein, erhielt: {n}")
        if n == 1:
            return 0

        aktuell = n
        schritte = 0

        # Schritt für Schritt bis zur 1
        while aktuell != 1:
            aktuell = self.step(aktuell)
            schritte += 1
            # Sicherheitsabbruch (Vermutung: immer endlich)
            if schritte > 10_000_000:
                raise RuntimeError(f"Maximale Schritte erreicht für n={n}")

        return schritte

    def trajectory_statistics(self, n: int) -> dict:
        """
        Berechnet Statistiken über die Collatz-Trajektorie von n.

        Gibt zurück:
        - max_value:          Maximaler Wert in der Folge
        - steps_to_1:         Anzahl Schritte bis 1 (totale Stoppzeit)
        - odd_steps:          Anzahl ungerader Schritte (3n+1 angewandt)
        - even_steps:         Anzahl gerader Schritte (n/2 angewandt)
        - compression_ratio:  steps_to_1 / log2(n) (relatives Maß)

        @param n: Startwert
        @return:  Dictionary mit Trajektorie-Statistiken
        """
        folge = self.sequence(n)

        ungerade = 0
        gerade = 0

        # Ungerade/gerade Schritte zählen (ohne letztes Element 1)
        for wert in folge[:-1]:
            if wert % 2 != 0:
                ungerade += 1
            else:
                gerade += 1

        schritte = len(folge) - 1  # Anzahl Übergänge
        log_n = math.log2(n) if n > 1 else 1.0

        return {
            "max_value":        max(folge),
            "steps_to_1":       schritte,
            "odd_steps":        ungerade,
            "even_steps":       gerade,
            "compression_ratio": schritte / log_n,
        }


# ===========================================================================
# Klasse: CollatzTwoadicAnalysis
# ===========================================================================

class CollatzTwoadicAnalysis:
    """
    Collatz als 2-adische Funktion (nach Lagarias).

    In der 2-adischen Sichtweise ist die Collatz-Funktion eine Kontraktion
    auf dem Raum der 2-adischen Zahlen ℤ₂. Der Lyapunov-Exponent

        λ = (1/2)·log(3/4) + (1/2)·log(1/2)  ≈ −0.0959

    ist negativ, was die beobachtete Konvergenz erklärt.

    Referenz: Lagarias (1985), "The 3x+1 problem and its generalizations"

    @author Michael Fuhrmann
    @since 2026-03-11
    """

    def __init__(self) -> None:
        """Initialisiert die 2-adische Collatz-Analyse."""
        self._collatz = CollatzFunction()

    def two_adic_valuation(self, n: int) -> int:
        """
        Berechnet die 2-adische Bewertung ν₂(n).

        ν₂(n) ist der größte Exponent k, so dass 2^k | n gilt.

        Beispiele:
            ν₂(8) = 3  (da 8 = 2³)
            ν₂(6) = 1  (da 6 = 2·3)
            ν₂(1) = 0  (da 1 = 2⁰·1)

        @param n: Positive ganze Zahl
        @return:  ν₂(n)
        @raises ValueError: Falls n = 0
        """
        if n == 0:
            raise ValueError("ν₂(0) ist undefiniert (formal +∞)")

        # Negativen Anteil ignorieren
        n = abs(n)
        k = 0
        # Solange n durch 2 teilbar: erhöhe Bewertung
        while n % 2 == 0:
            k += 1
            n //= 2
        return k

    def collatz_as_padic(self, n: int, steps: int) -> list:
        """
        Berechnet die Folge der 2-adischen Bewertungen ν₂(T^k(n)) für k=0..steps.

        Diese Folge zeigt, wie oft der aktuelle Wert durch 2 teilbar ist –
        ein Maß für die "2-adische Nähe zu 0".

        @param n:     Startwert
        @param steps: Anzahl Schritte
        @return:      Liste [ν₂(n), ν₂(T(n)), ν₂(T²(n)), ..., ν₂(T^steps(n))]
        """
        bewertungen = []
        aktuell = n

        for _ in range(steps + 1):
            bewertungen.append(self.two_adic_valuation(aktuell))
            if aktuell == 1:
                break
            aktuell = self._collatz.step(aktuell)

        return bewertungen

    def odd_part(self, n: int) -> int:
        """
        Berechnet den ungeraden Anteil von n.

        Der ungerade Anteil ist n / 2^{ν₂(n)}, also n mit allen
        Faktoren 2 herausgekürzt.

        Beispiele:
            odd_part(12) = 3  (12 = 4·3)
            odd_part(7)  = 7  (7 ungerade)

        @param n: Positive ganze Zahl
        @return:  Ungerader Anteil n / 2^{ν₂(n)}
        """
        n = abs(n)
        # Alle Zweierpotenzen herausteilen
        while n % 2 == 0:
            n //= 2
        return n

    def lyapunov_exponent(self, n: int, steps: int = 1000) -> float:
        """
        Schätzt den Lyapunov-Exponenten der Collatz-Folge für Startwert n.

        Definition:
            λ = lim_{k→∞} (1/k) · Σ_{j=0}^{k-1} log|T(T^j(n)) / T^j(n)|

        Für die Collatz-Funktion gilt:
            log(T(m)/m) = log(1/2)      falls m gerade
            log(T(m)/m) = log(3 + 1/m)  falls m ungerade  ≈ log(3)

        Erwarteter Wert (nach Lagarias):
            λ ≈ (1/2)·log(3) + (1/2)·log(1/2) = log(3)/2 − log(2) ≈ −0.0959

        Ein negativer Lyapunov-Exponent bedeutet: Abstände schrumpfen → Konvergenz.

        @param n:     Startwert
        @param steps: Anzahl Schritte für die Schätzung
        @return:      Geschätzter Lyapunov-Exponent λ
        """
        aktuell = n
        log_summe = 0.0
        gezaehlt = 0

        for _ in range(steps):
            if aktuell <= 1:
                break

            naechster = self._collatz.step(aktuell)
            # Logarithmischer Wachstumsfaktor dieses Schritts
            log_summe += math.log(naechster / aktuell)
            aktuell = naechster
            gezaehlt += 1

        if gezaehlt == 0:
            return 0.0

        return log_summe / gezaehlt

    def average_contraction(self, sample_size: int = 1000) -> float:
        """
        Schätzt die mittlere Kontraktion über vollständige Trajektorien.

        Berechnet für sample_size zufällige Startwerte den durchschnittlichen
        logarithmischen Wachstumsfaktor pro Schritt über die gesamte Trajektorie:

            λ̄ = (1/sample_size) · Σ_n lyapunov_exponent(n, steps=200)

        Ein negativer Wert bestätigt die Konvergenz der Collatz-Folge.
        Theoretisch: λ ≈ log(3)/2 − log(2) ≈ −0.0959

        @param sample_size: Anzahl zufällig gewählter Startwerte
        @return:            Empirische mittlere Kontraktion (erwartet: < 0)
        """
        log_summe = 0.0
        gezaehlt = 0

        # Zufällige Stichprobe von Startwerten, Lyapunov-Exponent je berechnen
        for _ in range(sample_size):
            n = random.randint(3, sample_size * 10)
            lam = self.lyapunov_exponent(n, steps=200)
            log_summe += lam
            gezaehlt += 1

        return log_summe / gezaehlt if gezaehlt > 0 else 0.0


# ===========================================================================
# Klasse: CollatzErgodicMeasure
# ===========================================================================

class CollatzErgodicMeasure:
    """
    Ergodische Maße auf dem Collatz-System (nach Terenz Tao, 2019).

    Taos Hauptresultat (2019): Für jede Funktion f: ℕ → ℝ mit f(n) → ∞ gilt:
        |{n ≤ N : min_{k≥0} T^k(n) ≤ f(n)}| / N → 1  (N → ∞)

    Das bedeutet: Fast alle Startwerte haben eine Trajektorie, die beliebig
    tief unter jede divergierende Schranke fällt.

    Referenz: Tao (2019), arXiv:1909.03562

    @author Michael Fuhrmann
    @since 2026-03-11
    """

    def __init__(self) -> None:
        """Initialisiert die ergodische Maß-Analyse."""
        self._collatz = CollatzFunction()

    def invariant_measure_estimate(self, N: int) -> np.ndarray:
        """
        Schätzt das empirische invariante Maß des Collatz-Systems.

        Berechnet die Häufigkeit, mit der jeder Wert in {1, ..., N} in
        allen Collatz-Folgen für n ≤ N vorkommt. Das Ergebnis ist normiert.

        Das theoretische invariante Maß konzentriert sich auf kleine Zahlen
        (Potenzen von 2, die zum Zyklus 1→4→2→1 gehören).

        @param N: Obere Grenze (berechnet Folgen für n = 1, ..., N)
        @return:  Normiertes empirisches Maß als numpy-Array der Länge N+1
                  (Index i enthält die relative Häufigkeit des Werts i)
        """
        # Häufigkeiten initialisieren (Index 0 unbenutzt)
        haeufigkeit = np.zeros(N + 1, dtype=float)

        for n in range(1, N + 1):
            aktuell = n
            # Collatz-Folge durchlaufen und Werte ≤ N zählen
            for _ in range(10000):
                if 1 <= aktuell <= N:
                    haeufigkeit[aktuell] += 1
                if aktuell == 1:
                    break
                aktuell = self._collatz.step(aktuell)

        # Normierung: Gesamthäufigkeit = 1
        gesamt = haeufigkeit.sum()
        if gesamt > 0:
            haeufigkeit /= gesamt

        return haeufigkeit

    def ergodic_average(self, f: Callable, N: int) -> float:
        """
        Berechnet den ergodischen Mittelwert einer messbaren Funktion f.

        Definition:
            Ā(f) = (1/N) · Σ_{n=1}^N f(T(n))

        Für ergodische Systeme konvergiert dies gegen den Raumittelwert
        ∫ f dμ bzgl. des invarianten Maßes μ.

        @param f: Messbare Funktion f: ℕ → ℝ
        @param N: Obere Grenze der Summe
        @return:  Ergodischer Mittelwert (1/N) Σ f(T(n))
        """
        summe = 0.0
        for n in range(1, N + 1):
            # Einen Collatz-Schritt von n aus
            naechster = self._collatz.step(n)
            summe += f(naechster)
        return summe / N

    def density_below_bound(self, N: int, bound_func: Callable) -> float:
        """
        Berechnet den Anteil der n ≤ N, deren Trajektorie die Schranke unterschreitet.

        Taos Theorem (grob): Für f(n) → ∞ gilt dieser Anteil → 1.

        @param N:          Obere Grenze
        @param bound_func: Schranken-Funktion f: ℕ → ℝ (z.B. n^ε)
        @return:           Anteil der n ≤ N mit min_k T^k(n) ≤ f(n)
        """
        treffer = 0

        for n in range(1, N + 1):
            schranke = bound_func(n)
            aktuell = n
            unterschritten = False

            # Trajektorie verfolgen
            for _ in range(10000):
                if aktuell <= schranke:
                    unterschritten = True
                    break
                if aktuell == 1:
                    # 1 ≤ schranke(n) für alle n ≥ 1 bei sinnvollen Schranken
                    unterschritten = True
                    break
                aktuell = self._collatz.step(aktuell)

            if unterschritten:
                treffer += 1

        return treffer / N

    def tao_theorem_numerical(self, N: int = 1000, epsilon: float = 0.5) -> dict:
        """
        Numerische Verifikation von Taos Satz für n ≤ N.

        Taos Satz: Für jedes ε > 0 gilt für fast alle n:
            min_{k≥0} T^k(n) ≤ n^ε

        Diese Funktion zählt, wie viele n ≤ N tatsächlich eine Trajektorie
        haben, die unter n^ε fällt.

        @param N:       Obere Grenze der Verifikation
        @param epsilon: Exponent der Schranke (typisch: 0 < ε < 1)
        @return:        Dictionary mit:
                        - "verified":  Anzahl bestätigter n
                        - "fraction":  Anteil bestätigter n (≈ 1 laut Tao)
                        - "epsilon":   Verwendeter Exponent ε
                        - "N":         Obere Grenze N
        """
        schranke_func = lambda n: n ** epsilon
        anteil = self.density_below_bound(N, schranke_func)
        anzahl = round(anteil * N)

        return {
            "verified": anzahl,
            "fraction": anteil,
            "epsilon":  epsilon,
            "N":        N,
        }


# ===========================================================================
# Klasse: CollatzStatistics
# ===========================================================================

class CollatzStatistics:
    """
    Statistische Analysen der Collatz-Folgen.

    Untersucht die Verteilung von Stoppzeiten, Rekordhalter und
    visualisiert Trajektorien und den Collatz-Baum.

    @author Michael Fuhrmann
    @since 2026-03-11
    """

    def __init__(self) -> None:
        """Initialisiert die Collatz-Statistik."""
        self._collatz = CollatzFunction()

    def stopping_time_distribution(self, limit: int) -> dict:
        """
        Berechnet die Verteilung der totalen Stoppzeiten für n ≤ limit.

        @param limit: Obere Grenze (berechnet Stoppzeiten für n = 1..limit)
        @return:      Dictionary {stoppzeit: anzahl_zahlen_mit_dieser_stoppzeit}
        """
        verteilung: dict[int, int] = {}

        for n in range(1, limit + 1):
            t = self._collatz.total_stopping_time(n)
            verteilung[t] = verteilung.get(t, 0) + 1

        return verteilung

    def record_setters(self, limit: int) -> list:
        """
        Findet alle Zahlen n ≤ limit mit rekordverdächtigen Stoppzeiten.

        Eine Zahl n setzt einen Stoppzeit-Rekord, falls ihre totale Stoppzeit
        größer ist als die aller kleineren Zahlen.

        @param limit: Obere Grenze
        @return:      Liste von (n, stoppzeit)-Tupeln in aufsteigender Reihenfolge
        """
        rekorde: list[tuple[int, int]] = []
        bisheriger_rekord = -1

        for n in range(1, limit + 1):
            t = self._collatz.total_stopping_time(n)
            if t > bisheriger_rekord:
                bisheriger_rekord = t
                rekorde.append((n, t))

        return rekorde

    def plot_stopping_times(self, limit: int = 1000) -> plt.Figure:
        """
        Erstellt einen Plot der totalen Stoppzeit τ(n) gegen n.

        @param limit: Obere Grenze der x-Achse
        @return:      matplotlib-Figure mit dem Plot
        """
        ns = list(range(1, limit + 1))
        stoppzeiten = [self._collatz.total_stopping_time(n) for n in ns]

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.scatter(ns, stoppzeiten, s=0.5, alpha=0.6, color="steelblue")
        ax.set_xlabel("n")
        ax.set_ylabel("Stoppzeit τ(n)")
        ax.set_title(f"Totale Stoppzeiten der Collatz-Folge für n ≤ {limit}")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        return fig

    def plot_trajectory(self, n: int) -> plt.Figure:
        """
        Erstellt einen logarithmischen Plot der Collatz-Trajektorie von n.

        @param n: Startwert
        @return:  matplotlib-Figure mit der Trajektorie (log-Skala y-Achse)
        """
        folge = self._collatz.sequence(n)
        schritte = list(range(len(folge)))

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(schritte, folge, color="darkorange", linewidth=0.8)
        ax.set_yscale("log")
        ax.set_xlabel("Schritt k")
        ax.set_ylabel("T^k(n) (log-Skala)")
        ax.set_title(f"Collatz-Trajektorie von n = {n} ({len(folge)-1} Schritte)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        return fig

    def collatz_tree(self, depth: int = 5) -> dict:
        """
        Erstellt den Collatz-Baum rückwärts: Welche Zahlen führen zu 1 in ≤ depth Schritten?

        Rückwärts von 1: Jede Zahl m hat Vorgänger:
            - 2m          (immer, da T(2m) = m)
            - (m-1)/3     falls (m-1) ≡ 0 mod 3 und (m-1)/3 ungerade

        @param depth: Tiefe des Baums (Schritte rückwärts von 1)
        @return:      Dictionary {zahl: vorgaenger_liste}
        """
        # Baum als Adjazenzliste (Eltern → Kinder im Rückwärtsgraph)
        baum: dict[int, list[int]] = {}
        # BFS rückwärts von 1
        aktuelle_ebene = {1}

        for _ in range(depth):
            naechste_ebene: set[int] = set()
            for m in aktuelle_ebene:
                vorgaenger: list[int] = []

                # Vorgänger 1: T(2m) = m → 2m ist immer Vorgänger
                v1 = 2 * m
                vorgaenger.append(v1)
                naechste_ebene.add(v1)

                # Vorgänger 2: T((m-1)/3) = m, falls (m-1) % 3 == 0
                # und (m-1)/3 ungerade (damit T(x) = 3x+1 gilt)
                if (m - 1) % 3 == 0:
                    v2 = (m - 1) // 3
                    if v2 > 0 and v2 % 2 != 0:  # Muss ungerade sein
                        vorgaenger.append(v2)
                        naechste_ebene.add(v2)

                baum[m] = vorgaenger

            aktuelle_ebene = naechste_ebene

        return baum


# ===========================================================================
# Freie Funktion: verify_collatz
# ===========================================================================

def verify_collatz(limit: int) -> tuple:
    """
    Verifiziert die Collatz-Vermutung für alle n ≤ limit.

    Prüft für jedes n = 1, 2, ..., limit, ob die Collatz-Folge die Zahl 1
    erreicht. Verwendet einen Cache für bereits verifizierte Zahlen.

    @param limit: Obere Grenze der Verifikation
    @return:      (True, limit)    falls alle n ≤ limit konvergieren
                  (False, n)       falls n das erste Gegenbeispiel ist
    """
    cf = CollatzFunction()

    # Cache: Zahlen die bereits als konvergent bekannt sind
    # (1 ist trivial konvergent)
    konvergent: set[int] = {1}

    for n in range(1, limit + 1):
        if n in konvergent:
            continue

        # Trajektorie verfolgen
        trajektorie = []
        aktuell = n

        while aktuell not in konvergent:
            if aktuell == 1:
                break
            trajektorie.append(aktuell)
            aktuell = cf.step(aktuell)

            # Sicherheitsabbruch: unerwarteter langer Orbit
            if len(trajektorie) > 10_000_000:
                return (False, n)

        # Alle Werte in der Trajektorie als konvergent markieren
        for wert in trajektorie:
            konvergent.add(wert)

    return (True, limit)
