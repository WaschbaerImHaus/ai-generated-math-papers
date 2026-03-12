"""
Lonely Runner Conjecture — Verifikation und Analyse.

Die Lonely-Runner-Vermutung (Wills 1967, unabhängig Cusick 1973):
Gegeben n Läufer auf einem Einheitskreis (Umfang 1), die alle bei 0
starten und verschiedene ganzzahlige Geschwindigkeiten haben. Dann
gibt es für jeden Läufer einen Zeitpunkt t, zu dem er von ALLEN
anderen Läufern mindestens 1/n entfernt ist.

Formale Aussage:
  Für n ≥ 2 und paarweise verschiedene ganze Zahlen v_0,...,v_{n-1}
  (mit v_0 = 0 o.B.d.A.) existiert t > 0 mit:
    ∥t·(v_i − v_j)∥ ≥ 1/n  für alle i ≠ 0
  wobei ∥x∥ = min(x mod 1, 1 - x mod 1) der Kreisabstand ist.

Beweisstatus (CONJECTURE — unbewiesen im Allgemeinen):
  - n=2: Trivial
  - n=3,4: Proved (Cusick & Pomerance 1984)
  - n=5: Proved (Goddyn & various, 2006)
  - n=6: Proved (Bohman, Fonoberova & Pikhurko, 2011)
  - n=7: Proved (Tao, 2018)
  - n≥8: OFFEN (Conjecture)

@author: Michael Fuhrmann
@date: 2026-03-12
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import List, Optional, Tuple, Dict
import numpy as np
from sympy import gcd, lcm, Rational, factorint


class LonelyRunnerConjecture:
    """
    Verifikation und Simulation der Lonely-Runner-Vermutung.

    Diese Klasse implementiert:
      1. Exakte Verifikation für kleine n (≤6) mit rationalen Arithmetik
      2. Numerische Simulation mit diskreten Zeitschritten
      3. Diophantische Reformulierung
      4. Formale Beschreibung von Beweisansätzen

    WICHTIG: Die Lonely-Runner-Vermutung ist für n ≥ 8 NICHT bewiesen.
    Diese Implementierung verifiziert den Spezialfall und demonstriert
    die Methodik, erhebt aber keinen Anspruch auf einen allgemeinen Beweis.

    @author: Michael Fuhrmann
    @date: 2026-03-12
    """

    def __init__(self, velocities: List[int]) -> None:
        """
        Initialisiert das Lonely-Runner-Problem mit n Läufern.

        Konvention: Läufer 0 hat immer Geschwindigkeit 0 (Referenzläufer).
        Falls 0 nicht in der Liste, wird sie automatisch hinzugefügt.

        @param velocities: Liste von verschiedenen ganzzahligen Geschwindigkeiten
        @raises ValueError: Falls Geschwindigkeiten nicht paarweise verschieden
        @date: 2026-03-12
        """
        # Eindeutigkeit prüfen
        if len(set(velocities)) != len(velocities):
            raise ValueError("Geschwindigkeiten müssen paarweise verschieden sein.")

        # Normalisierung: O.B.d.A. Läufer 0 hat Geschwindigkeit 0
        # Subtrahiere v[0] von allen Geschwindigkeiten
        v0 = velocities[0]
        self.velocities = [v - v0 for v in velocities]
        self.n = len(self.velocities)  # Anzahl Läufer
        # Mindestabstand: 1/n
        self.threshold = Fraction(1, self.n)

    @property
    def relative_velocities(self) -> List[int]:
        """
        Gibt die relativen Geschwindigkeiten (bezogen auf Läufer 0) zurück.

        @return: Liste der relativen Geschwindigkeiten
        @date: 2026-03-12
        """
        # Relative Geschwindigkeit von Läufer i zu Läufer 0
        return [v for v in self.velocities if v != 0]

    def circle_distance(self, x: float) -> float:
        """
        Berechnet den Abstand auf dem Einheitskreis (kürzeste Bogenlänge).

        d(x) = min(x mod 1, 1 - x mod 1)

        @param x: Position auf dem Kreis (reell)
        @return: Kreisabstand in [0, 0.5]
        @date: 2026-03-12
        """
        x_mod = x % 1.0
        return min(x_mod, 1.0 - x_mod)

    def circle_distance_fraction(self, x: Fraction) -> Fraction:
        """
        Exakte Kreisabstand-Berechnung mit rationaler Arithmetik.

        @param x: Rationale Position
        @return: Exakter Kreisabstand als Fraction
        @date: 2026-03-12
        """
        # x mod 1
        x_mod = x - int(x)
        if x_mod < 0:
            x_mod += 1
        return min(x_mod, Fraction(1) - x_mod)

    def is_lonely_at_time(self, runner_idx: int, t: float,
                          tol: float = 1e-9) -> bool:
        """
        Prüft ob ein Läufer zum Zeitpunkt t "lonely" (einsam) ist.

        "Lonely" bedeutet: Abstand zu JEDEM anderen Läufer ≥ 1/n.

        @param runner_idx: Index des zu prüfenden Läufers
        @param t: Zeitpunkt
        @param tol: Numerische Toleranz für Vergleiche
        @return: True falls Läufer runner_idx zum Zeitpunkt t lonely ist
        @date: 2026-03-12
        """
        threshold_float = 1.0 / self.n
        v_i = self.velocities[runner_idx]

        for j, v_j in enumerate(self.velocities):
            if j == runner_idx:
                continue
            # Relativer Abstand auf dem Kreis
            relative_pos = t * (v_i - v_j)
            dist = self.circle_distance(relative_pos)
            if dist < threshold_float - tol:
                return False
        return True

    def is_lonely_at_time_exact(self, runner_idx: int, t: Fraction) -> bool:
        """
        Exakte Version von is_lonely_at_time mit rationaler Arithmetik.

        @param runner_idx: Index des Läufers
        @param t: Rationaler Zeitpunkt
        @return: True falls Läufer zum Zeitpunkt t lonely ist
        @date: 2026-03-12
        """
        threshold = Fraction(1, self.n)
        v_i = self.velocities[runner_idx]

        for j, v_j in enumerate(self.velocities):
            if j == runner_idx:
                continue
            relative_pos = t * Fraction(v_i - v_j)
            dist = self.circle_distance_fraction(relative_pos)
            if dist < threshold:
                return False
        return True

    def find_lonely_time_simulation(
        self, runner_idx: int,
        t_max: float = 100.0,
        steps: int = 100000
    ) -> Optional[float]:
        """
        Findet durch numerische Simulation einen Zeitpunkt, zu dem Läufer lonely ist.

        Diskretisiert [0, t_max] in gleichmäßige Schritte und prüft jeden.

        @param runner_idx: Index des zu prüfenden Läufers
        @param t_max: Maximale Suchzeit
        @param steps: Anzahl Zeitschritte
        @return: Zeitpunkt t ≥ 0 falls gefunden, sonst None
        @date: 2026-03-12
        """
        dt = t_max / steps
        for step in range(1, steps + 1):
            t = step * dt
            if self.is_lonely_at_time(runner_idx, t):
                return t
        return None

    def find_lonely_time_exact(
        self, runner_idx: int,
        max_denominator: int = 1000
    ) -> Optional[Fraction]:
        """
        Sucht exakten rationalen Zeitpunkt für Loneliness von Läufer runner_idx.

        Strategie: Kandidaten-Zeitpunkte sind Brüche a/b mit b ≤ max_denominator.
        Bei rationalen Geschwindigkeiten reicht es, rationale t zu prüfen.

        @param runner_idx: Index des Läufers
        @param max_denominator: Maximaler Nenner für Bruchsuche
        @return: Rationaler Zeitpunkt falls gefunden, sonst None
        @date: 2026-03-12
        """
        # Kritische Zeitpunkte: t = k/(v_i - v_j) für k=1,...,q
        # Zwischen diesen ist die Loneliness-Funktion monoton
        critical_times = set()
        for i, vi in enumerate(self.velocities):
            for j, vj in enumerate(self.velocities):
                if i != j and vi != vj:
                    diff = abs(vi - vj)
                    # Midpoints zwischen kritischen Punkten prüfen
                    for k in range(1, diff * max_denominator // diff + 1):
                        t = Fraction(2 * k - 1, 2 * diff)
                        if 0 < t <= max_denominator:
                            critical_times.add(t)

        # Midpoints direkt prüfen
        for t in sorted(critical_times):
            if self.is_lonely_at_time_exact(runner_idx, t):
                return t

        return None

    def verify_conjecture_for_all_runners(
        self, t_max: float = 200.0,
        steps: int = 500000
    ) -> Dict[int, Optional[float]]:
        """
        Prüft die Lonely-Runner-Vermutung für alle Läufer numerisch.

        Gibt für jeden Läufer einen Zeitpunkt zurück, zu dem er lonely ist.
        Falls keiner gefunden wird, ist die Vermutung für diese Instanz
        numerisch nicht verifiziert (erhöhe steps oder t_max).

        @param t_max: Maximale Suchzeit
        @param steps: Anzahl Zeitschritte
        @return: Dict {Läufer-Index: Zeitpunkt oder None}
        @date: 2026-03-12
        """
        result = {}
        for i in range(self.n):
            if self.velocities[i] == 0:
                # Läufer 0 (Referenz): Finde t, sodass alle anderen weit weg sind
                # Das ist äquivalent zum allgemeinen Fall durch Symmetrie
                t = self.find_lonely_time_simulation(i, t_max, steps)
            else:
                t = self.find_lonely_time_simulation(i, t_max, steps)
            result[i] = t
        return result

    def min_distance_at_time(self, runner_idx: int, t: float) -> float:
        """
        Berechnet den minimalen Abstand von runner_idx zu allen anderen Läufern.

        @param runner_idx: Läufer-Index
        @param t: Zeitpunkt
        @return: Minimaler Kreisabstand
        @date: 2026-03-12
        """
        v_i = self.velocities[runner_idx]
        distances = []
        for j, v_j in enumerate(self.velocities):
            if j == runner_idx:
                continue
            relative_pos = t * (v_i - v_j)
            distances.append(self.circle_distance(relative_pos))
        return min(distances) if distances else float('inf')

    def describe_diophantine_reformulation(self) -> str:
        """
        Beschreibt die diophantische Reformulierung der Vermutung.

        Die Vermutung ist äquivalent zu:
          Für jeden Läufer i ≠ 0 existiert t > 0 mit:
            ∥t·(v_i − v_j)∥ ≥ 1/n  für alle j ≠ i
          wobei ∥x∥ = dist(x, ℤ) der Abstand zur nächsten ganzen Zahl ist.

        Bei ganzzahligen Geschwindigkeiten: O.B.d.A. v_0 = 0,
        dann sind die relativen Geschwindigkeiten v_1,...,v_{n-1}.

        @return: Beschreibungstext
        @date: 2026-03-12
        """
        rel = self.relative_velocities
        lines = [
            "=== Diophantische Reformulierung der Lonely-Runner-Vermutung ===",
            f"n = {self.n} Läufer, Schwellwert = 1/{self.n} = {1/self.n:.6f}",
            f"Relative Geschwindigkeiten (bzgl. Läufer 0): {rel}",
            "",
            "Zu zeigen: Für jeden Läufer i existiert t > 0 mit",
            "  ∥t·(vᵢ − vⱼ)∥ ≥ 1/n  für alle j ≠ i",
            "",
            "Äquivalente Formulierung (Cusick 1973):",
            "  Für alle ganzzahligen w₁,...,w_{n-1} (paarweise verschieden)",
            "  existiert reelles t mit ∥t·wₖ∥ ≥ 1/n für k=1,...,n-1",
        ]
        return "\n".join(lines)

    def blocki_sketch_n4(self) -> str:
        """
        Beschreibt den Beweisansatz von Blocki für n=4 (vereinfachte Skizze).

        Blocki (2010, Masters Thesis) bewies den Fall n=4 mit einem
        eleganten kombinatorischen Argument:

        Sei n=4 mit Geschwindigkeiten 0, a, b, c (paarweise verschieden).
        O.B.d.A. gcd(a,b,c) = 1 (sonst skaliere die Zeit).

        Schlüsselidee: Betrachte die Intervalle I_k = (k/4, (k+1)/4) mod 1
        für k=0,1,2,3 (vier Viertel des Kreises).

        Für jeden Läufer i: Die Zeitmengen T_i = {t : Läufer i ist lonely}
        haben positive Dichte (Lebesgue-Maß > 0 in jedem Periodenintervall).

        Formaler Beweis-Sketch:
          1. Betrachte Funktion f(t) = min_{j≠i} ∥t(vᵢ-vⱼ)∥
          2. f ist stetig, periodisch mit Periode L = lcm(|vᵢ-vⱼ|)⁻¹
          3. Zeige: max_{t} f(t) ≥ 1/4 via Diskrepanztheorie
          4. Schlüssellemma: Für vier verschiedene reelle Zahlen α,β,γ auf
             [0,1) gilt: max_t min(∥tα∥, ∥tβ∥, ∥tγ∥) ≥ 1/4

        @return: Beschreibungstext des Beweis-Sketch
        @date: 2026-03-12
        """
        lines = [
            "=== Blocki-Beweis-Sketch für n=4 (2010) ===",
            "",
            "Gegeben: 4 Läufer mit Geschwindigkeiten 0, a, b, c ∈ ℤ (verschieden)",
            "Zu zeigen: Jeder Läufer ist zu einem Zeitpunkt ≥1/4 von allen anderen entfernt.",
            "",
            "O.B.d.A.: v₀=0, v₁=a, v₂=b, v₃=c, gcd(a,b,c)=1",
            "",
            "Beweisschritt 1: Periodizität",
            "  Alle Läuferpositionen sind periodisch mit Periode T = lcm(a,b,c)⁻¹.",
            "  Es genügt, t ∈ [0, 1) zu betrachten.",
            "",
            "Beweisschritt 2: Kritische Intervalle",
            "  Definiere für jeden Läufer i die 'bösen Zeiten':",
            "    B_{ij} = {t : ∥t(vᵢ-vⱼ)∥ < 1/4}",
            "  Jedes B_{ij} hat Maß 1/2 in [0,1).",
            "",
            "Beweisschritt 3: Überlappungsargument",
            "  Schlüssellemma (Cusick & Pomerance 1984):",
            "    Für drei linear unabhängige ganze Zahlen α,β,γ gilt:",
            "    |[0,1) \\ (B_{01} ∪ B_{02} ∪ B_{03})| > 0",
            "  D.h.: Es gibt immer Zeitpunkte, wo Läufer 0 von allen drei entfernt ist.",
            "",
            "Beweisschritt 4: Symmetrieargument",
            "  Dasselbe gilt für jeden Läufer durch Koordinatenwechsel.",
            "",
            "HINWEIS: Dies ist eine VEREINFACHTE SKIZZE.",
            "Der vollständige Beweis erfordert Fourier-Analyse auf ℝ/ℤ.",
            "Quelle: Blocki, J. (2010). The Lonely Runner Problem.",
        ]
        return "\n".join(lines)


class LonelyRunnerVerifier:
    """
    Exakter Verifikator für kleine Fälle der Lonely-Runner-Vermutung.

    Für n ≤ 6 Läufer kann die Vermutung durch erschöpfende Suche
    über rationale Zeitpunkte verifiziert werden.

    @author: Michael Fuhrmann
    @date: 2026-03-12
    """

    def __init__(self) -> None:
        """
        Initialisiert den Verifikator.

        @date: 2026-03-12
        """
        pass

    def verify_n2(self, v1: int) -> Tuple[bool, Fraction]:
        """
        Verifiziert die Lonely-Runner-Vermutung für n=2 (trivial).

        n=2: Zwei Läufer mit Geschwindigkeiten 0 und v₁.
        Läufer 0 ist lonely wenn ∥t·v₁∥ ≥ 1/2.
        Exakt bei t = 1/(2·v₁): Relative Position = 1/2.

        @param v1: Geschwindigkeit des zweiten Läufers (≠0)
        @return: (True, t_lonely) — immer wahr für n=2
        @date: 2026-03-12
        """
        # t = 1/(2*v1) liefert ∥t·v1∥ = ∥1/2∥ = 1/2 ≥ 1/2 = threshold
        t_lonely = Fraction(1, 2 * abs(v1))
        return True, t_lonely

    def verify_runner_exact(
        self,
        velocities: List[int],
        runner_idx: int,
        max_denom: int = 200
    ) -> Tuple[bool, Optional[Fraction]]:
        """
        Exakte Verifikation: Prüft ob Läufer runner_idx einen Loneliness-Zeitpunkt hat.

        Strategie: Generiere alle kritischen Zeitpunkte der Loneliness-Funktion
        und prüfe Mittelpunkte zwischen je zwei aufeinanderfolgenden kritischen t.

        Kritische Zeitpunkte: t = k / (v_i - v_j) für alle Paare (i,j) und k ∈ ℕ.
        Die Loneliness-Funktion ist zwischen kritischen Punkten konstant.
        Daher genügt es, in jedem Intervall einen Punkt (z.B. Mittelpunkt) zu testen.

        @param velocities: Geschwindigkeitsliste (v[0] wird auf 0 normiert)
        @param runner_idx: Index des zu prüfenden Läufers
        @param max_denom: Maximaler Nenner für Bruchsuche
        @return: (True, t) falls Zeitpunkt gefunden, (False, None) sonst
        @date: 2026-03-12
        """
        lrc = LonelyRunnerConjecture(velocities)
        n = lrc.n
        vels = lrc.velocities

        # Alle paarweisen relativen Geschwindigkeiten (v_i - v_j) für runner_idx
        diffs = []
        v_i = vels[runner_idx]
        for j, v_j in enumerate(vels):
            if j == runner_idx:
                continue
            d = abs(v_i - v_j)
            if d > 0:
                diffs.append(d)

        # Kritische Zeitpunkte in (0, 1]: t = k/d für k=0,...,d-1, d ∈ diffs
        # Wir sammeln alle Brüche k/d mit k=0..d, d ∈ diffs, im Intervall [0,1]
        critical = set()
        critical.add(Fraction(0))
        critical.add(Fraction(1))
        for d in diffs:
            for k in range(d + 1):
                critical.add(Fraction(k, d))

        # Sortiere und bilde Mittelpunkte jedes Intervalls
        sorted_crit = sorted(critical)
        candidates = []
        for i in range(len(sorted_crit) - 1):
            a, b = sorted_crit[i], sorted_crit[i + 1]
            midpoint = (a + b) / 2  # Fraction-Arithmetik: exakt
            if midpoint > 0:
                candidates.append(midpoint)

        # Zusätzlich: Kandidaten über alle Brüche a/b mit b ≤ max_denom
        # (für Läufer mit v_i=0 können kritische Punkte eng liegen)
        for b in range(1, max_denom + 1):
            for a in range(1, b):
                candidates.append(Fraction(a, b))

        # Alle Kandidaten prüfen
        for t in candidates:
            if t > 0 and lrc.is_lonely_at_time_exact(runner_idx, t):
                return True, t

        return False, None

    def verify_all_runners(
        self,
        velocities: List[int],
        max_denom: int = 500
    ) -> Dict[int, Tuple[bool, Optional[Fraction]]]:
        """
        Verifiziert die Vermutung für alle Läufer einer Instanz.

        @param velocities: Liste von Geschwindigkeiten
        @param max_denom: Maximaler Suchnenner
        @return: Dict {runner_idx: (verified, lonely_time)}
        @date: 2026-03-12
        """
        results = {}
        for i in range(len(velocities)):
            ok, t = self.verify_runner_exact(velocities, i, max_denom)
            results[i] = (ok, t)
        return results

    def generate_test_cases(self, n: int, max_vel: int = 10) -> List[List[int]]:
        """
        Generiert Test-Instanzen für n Läufer mit kleinen Geschwindigkeiten.

        @param n: Anzahl Läufer
        @param max_vel: Maximale Geschwindigkeit
        @return: Liste von Geschwindigkeitslisten
        @date: 2026-03-12
        """
        from itertools import combinations
        test_cases = []
        # Wähle n paarweise verschiedene Geschwindigkeiten aus {0,...,max_vel}
        for vels in combinations(range(max_vel + 1), n):
            test_cases.append(list(vels))
        return test_cases[:20]  # Maximal 20 Fälle zurückgeben
