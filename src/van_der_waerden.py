"""
Van der Waerden Zahlen und AP-Färbungsprobleme.

Dieses Modul implementiert die Berechnung und Verifikation von Van-der-Waerden-Zahlen
W(k;r): Die kleinste Zahl N, sodass jede r-Färbung von {1,...,N} eine
monochromatische arithmetische Progression (AP) der Länge k enthält.

Bekannte exakte Werte:
  W(2;3) = 9,  W(2;4) = 35,  W(2;5) = 178
  W(3;3) = 27, W(3;4) = 293

Mathematische Grundlage:
  Van der Waerdens Theorem (1927): Für alle k,r ≥ 1 existiert W(k;r).
  Gowers (2001) liefert obere Schranken über Fourier-Methoden:
    W(k;2) ≤ 2^{2^{2^{...}}} (Türme von 2, Primzahl-freie Variante)
  Berlekamp (1968): W(p+1;2) > p·2^p für jede Primzahl p.

@author: Michael Fuhrmann
@date: 2026-03-12
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from typing import List, Optional, Tuple, Dict, Set

import numpy as np
from pysat.solvers import Solver
from pysat.formula import CNF


class ArithmeticProgressionColoring:
    """
    SAT-Kodierung des AP-Färbungsproblems.

    Kodiert die Frage: "Gibt es eine r-Färbung von {1,...,N},
    die KEINE monochromatische k-term AP enthält?" als CNF-Formel.

    SAT-Variablen: x_{i,c} = True bedeutet, Zahl i erhält Farbe c.
    (interne Variable-ID: (i-1)*r + c,  1-basiert)

    @author: Michael Fuhrmann
    @date: 2026-03-12
    """

    def __init__(self, N: int, k: int, r: int) -> None:
        """
        Initialisiert die SAT-Kodierung.

        @param N: Größe des Bereichs {1,...,N}
        @param k: Länge der gesuchten AP
        @param r: Anzahl Farben
        @date: 2026-03-12
        """
        # Eingabeparameter speichern
        self.N = N
        self.k = k
        self.r = r
        # Alle APs vorberechnen (effizienter als on-the-fly)
        self._aps: List[List[int]] = self._compute_all_aps()

    def _var(self, i: int, c: int) -> int:
        """
        Berechnet die SAT-Variablen-ID für (Zahl i, Farbe c).

        Variablen sind 1-basiert (pysat-Konvention).
        Formel: var(i,c) = (i-1)*r + c

        @param i: Zahl aus {1,...,N}
        @param c: Farbe aus {1,...,r}
        @return: SAT-Variablen-ID (≥1)
        @date: 2026-03-12
        """
        return (i - 1) * self.r + c

    def _compute_all_aps(self) -> List[List[int]]:
        """
        Berechnet alle arithmetischen Progressionen der Länge k in {1,...,N}.

        Eine AP der Länge k mit Startwert a und Schritt d ist:
          a, a+d, a+2d, ..., a+(k-1)d

        @return: Liste aller k-term APs in {1,...,N}
        @date: 2026-03-12
        """
        aps = []
        for a in range(1, self.N + 1):
            for d in range(1, self.N):
                # AP: a, a+d, ..., a+(k-1)*d muss ≤ N bleiben
                ap = [a + j * d for j in range(self.k)]
                if ap[-1] > self.N:
                    break  # größere d bringen nur größere Endwerte
                aps.append(ap)
        return aps

    def build_cnf(self) -> CNF:
        """
        Erzeugt die CNF-Formel für das Färbungsproblem.

        Drei Typen von Klauseln:
          1) Mindestens-eine-Farbe: ∨_c x_{i,c} für alle i
          2) Höchstens-eine-Farbe: ¬x_{i,c1} ∨ ¬x_{i,c2} für c1≠c2
          3) Keine monochromatische AP: ¬x_{a,c} ∨ ¬x_{a+d,c} ∨ ...
             (für alle APs und alle Farben c)

        @return: CNF-Formel als pysat.formula.CNF
        @date: 2026-03-12
        """
        cnf = CNF()

        # --- Typ 1: Jede Zahl bekommt mindestens eine Farbe ---
        for i in range(1, self.N + 1):
            cnf.append([self._var(i, c) for c in range(1, self.r + 1)])

        # --- Typ 2: Jede Zahl bekommt höchstens eine Farbe ---
        for i in range(1, self.N + 1):
            for c1 in range(1, self.r + 1):
                for c2 in range(c1 + 1, self.r + 1):
                    cnf.append([-self._var(i, c1), -self._var(i, c2)])

        # --- Typ 3: Keine monochromatische AP ---
        for ap in self._aps:
            for c in range(1, self.r + 1):
                # Für jede AP und jede Farbe: nicht alle Elemente gleiche Farbe
                clause = [-self._var(x, c) for x in ap]
                cnf.append(clause)

        return cnf

    def solve(self) -> Optional[List[int]]:
        """
        Löst das SAT-Problem: Gibt es eine gültige Färbung ohne monochromatische AP?

        @return: Färbung als Liste (Index i-1 → Farbe 1..r) falls ERFÜLLBAR,
                 None falls UNERFÜLLBAR (d.h. jede Färbung enthält eine AP)
        @date: 2026-03-12
        """
        cnf = self.build_cnf()
        with Solver(bootstrap_with=cnf) as solver:
            if not solver.solve():
                return None
            model = solver.get_model()

        # Färbung aus Modell extrahieren
        coloring = [0] * (self.N + 1)  # Index 0 ungenutzt, 1..N genutzt
        for i in range(1, self.N + 1):
            for c in range(1, self.r + 1):
                var_id = self._var(i, c)
                if var_id <= len(model) and model[var_id - 1] > 0:
                    coloring[i] = c
                    break
        return coloring[1:]  # Nur Einträge 1..N zurückgeben

    def extract_coloring_dict(self) -> Optional[Dict[int, int]]:
        """
        Gibt die Färbung als Dictionary {Zahl: Farbe} zurück.

        @return: Dict oder None falls keine gültige Färbung existiert
        @date: 2026-03-12
        """
        result = self.solve()
        if result is None:
            return None
        return {i + 1: result[i] for i in range(len(result))}


class VanDerWaerdenNumbers:
    """
    Berechnung und Verifikation von Van-der-Waerden-Zahlen W(k;r).

    W(k;r) ist die kleinste natürliche Zahl N, sodass jede r-Färbung
    von {1,...,N} eine monochromatische arithmetische Progression der
    Länge k enthält.

    Bekannte Werte (Tabelle aus Literatur):
      W(2;3) = 9
      W(2;4) = 35
      W(2;5) = 178
      W(3;3) = 27
      W(3;4) = 293

    Für große Parameter sind exakte Werte unbekannt.

    @author: Michael Fuhrmann
    @date: 2026-03-12
    """

    # Tabelle bekannter exakter Werte: (k, r) → W(k;r)
    # Konvention: k = AP-Länge, r = Anzahl Farben
    # W(3;2)=9: Kleinste N sodass jede 2-Färbung eine monochromatische 3-AP enthält
    KNOWN_VALUES: Dict[Tuple[int, int], int] = {
        (3, 2): 9,
        (4, 2): 35,
        (5, 2): 178,
        (3, 3): 27,
        (3, 4): 293,
    }

    def __init__(self) -> None:
        """
        Initialisiert den Van-der-Waerden-Zahlen-Rechner.

        @date: 2026-03-12
        """
        # Cache für berechnete Werte (verhindert Neuberechnung)
        self._cache: Dict[Tuple[int, int], int] = {}

    def lookup(self, k: int, r: int) -> Optional[int]:
        """
        Gibt den bekannten exakten Wert W(k;r) aus der Tabelle zurück.

        @param k: AP-Länge
        @param r: Anzahl Farben
        @return: W(k;r) falls bekannt, sonst None
        @date: 2026-03-12
        """
        return self.KNOWN_VALUES.get((k, r))

    def is_ap_free_coloring(self, coloring: List[int], k: int) -> bool:
        """
        Prüft ob eine gegebene Färbung keine monochromatische k-term AP enthält.

        @param coloring: Liste von Farben (1-basiert), Länge N
        @param k: Länge der gesuchten AP
        @return: True falls keine monochromatische AP vorhanden
        @date: 2026-03-12
        """
        N = len(coloring)
        for a in range(N):
            for d in range(1, N):
                # AP: a, a+d, ..., a+(k-1)d (0-basierte Indizes)
                indices = [a + j * d for j in range(k)]
                if indices[-1] >= N:
                    break
                # Prüfe ob alle gleiche Farbe haben
                colors = [coloring[idx] for idx in indices]
                if len(set(colors)) == 1:
                    return False
        return True

    def has_monochromatic_ap(self, coloring: List[int], k: int) -> bool:
        """
        Prüft ob eine Färbung eine monochromatische k-term AP enthält.

        @param coloring: Färbung als Liste (0-basiert)
        @param k: AP-Länge
        @return: True falls monochromatische AP vorhanden
        @date: 2026-03-12
        """
        return not self.is_ap_free_coloring(coloring, k)

    def find_ap_free_coloring(self, N: int, k: int, r: int) -> Optional[List[int]]:
        """
        Versucht eine r-Färbung von {1,...,N} ohne k-term AP zu finden.

        Verwendet den SAT-Solver für effiziente Suche.

        @param N: Bereichsgröße
        @param k: AP-Länge
        @param r: Anzahl Farben
        @return: Gültige Färbung (Liste, 0-basiert) oder None
        @date: 2026-03-12
        """
        coder = ArithmeticProgressionColoring(N, k, r)
        return coder.solve()

    def verify_lower_bound(self, W_value: int, k: int, r: int) -> bool:
        """
        Verifiziert die Untergrenze: Es gibt eine (W-1)-Färbung ohne AP.

        W(k;r) ≥ W_value ⟺ {1,...,W_value-1} ist r-färbbar ohne k-term AP.

        @param W_value: Zu verifizierende Untergrenze
        @param k: AP-Länge
        @param r: Anzahl Farben
        @return: True falls Untergrenze verifiziert
        @date: 2026-03-12
        """
        coloring = self.find_ap_free_coloring(W_value - 1, k, r)
        return coloring is not None

    def verify_upper_bound(self, W_value: int, k: int, r: int) -> bool:
        """
        Verifiziert die Obergrenze: Jede W_value-Färbung enthält eine AP.

        W(k;r) ≤ W_value ⟺ {1,...,W_value} ist NICHT r-färbbar ohne k-term AP.

        @param W_value: Zu verifizierende Obergrenze
        @param k: AP-Länge
        @param r: Anzahl Farben
        @return: True falls Obergrenze verifiziert
        @date: 2026-03-12
        """
        coloring = self.find_ap_free_coloring(W_value, k, r)
        return coloring is None

    def compute(self, k: int, r: int, max_N: int = 50) -> Optional[int]:
        """
        Berechnet W(k;r) durch schrittweise SAT-Suche.

        Algorithmus:
          Für N = k, k+1, k+2, ...:
            Falls SAT(N, k, r) = UNSAT: W(k;r) = N, fertig.
            Falls SAT(N, k, r) = SAT:   weiter suchen.

        Warnung: Läuft bis max_N, falls kein Ergebnis gefunden.
        Für große W(k;r) empfiehlt sich lookup().

        @param k: AP-Länge (≥2)
        @param r: Anzahl Farben (≥2)
        @param max_N: Maximaler Suchbereich
        @return: W(k;r) oder None falls nicht in [k, max_N] gefunden
        @date: 2026-03-12
        """
        # Cache-Lookup
        key = (k, r)
        if key in self._cache:
            cached = self._cache[key]
            # Respektiere max_N: Bekannter Wert > max_N → None (wie SAT-Suche)
            return cached if cached <= max_N else None
        # Bekannte Werte nur verwenden wenn sie im Suchbereich liegen
        if key in self.KNOWN_VALUES:
            known = self.KNOWN_VALUES[key]
            if known <= max_N:
                return known
            # Wert bekannt aber außerhalb max_N → trotzdem via SAT prüfen

        # Schrittweise SAT-Suche
        for N in range(k, max_N + 1):
            coloring = self.find_ap_free_coloring(N, k, r)
            if coloring is None:
                # Kein AP-freies Zeugen für N → W(k;r) = N
                self._cache[key] = N
                return N

        return None  # Innerhalb max_N nicht gefunden

    def berlekamp_lower_bound(self, p: int) -> int:
        """
        Berechnet die Berlekamp-Untergrenze: W(p+1; 2) > p·2^p.

        Berlekamp (1968): Für jede Primzahl p gilt W(p+1; 2) > p·2^p.
        Dies ergibt sich aus der expliziten Konstruktion einer AP-freien
        2-Färbung von {1,...,p·2^p} mittels Galois-Körper GF(2^p).

        @param p: Primzahl
        @return: Berlekamp-Untergrenze p·2^p
        @date: 2026-03-12
        """
        # Berlekamp-Untergrenze: p * 2^p
        return p * (2 ** p)

    def gowers_upper_bound(self, k: int) -> int:
        """
        Berechnet eine grobe Gowers-Schranke für W(k;2).

        Gowers (2001) zeigte W(k;2) ≤ 2^{2^{k+9}}.
        Diese Schranke ist astronomisch groß und hier nur zur
        Vollständigkeit implementiert.

        HINWEIS: Für k≥4 liefert diese Formel unhandliche Zahlen.
        Nur für kleine k nutzbar.

        @param k: AP-Länge
        @return: Gowers-Obergrenze 2^{2^{k+9}}
        @date: 2026-03-12
        """
        # Gowers-Schranke: 2^(2^(k+9))
        # Sehr schnell wachsend — Vorsicht bei k≥4
        exponent = 2 ** (k + 9)
        # Rückgabe als Integer (kann sehr groß werden)
        return 2 ** exponent

    def detect_ap(self, sequence: List[int], k: int) -> List[List[int]]:
        """
        Findet alle arithmetischen Progressionen der Länge k in einer Folge.

        Eine AP (a, a+d, a+2d, ...) muss als Teilfolge (nicht notwendig
        zusammenhängend) in sequence vorkommen.

        @param sequence: Zu durchsuchende Folge
        @param k: Gesuchte AP-Länge
        @return: Liste aller gefundenen APs
        @date: 2026-03-12
        """
        found = []
        n = len(sequence)
        for i in range(n):
            for j in range(i + 1, n):
                # Berechne gemeinsame Differenz d aus sequence[i] und sequence[j]
                d = sequence[j] - sequence[i]
                if d <= 0:
                    continue
                # Prüfe ob weitere k-2 Elemente folgen
                ap = [sequence[i] + m * d for m in range(k)]
                if all(x in sequence for x in ap):
                    if ap not in found:
                        found.append(ap)
        return found

    def is_arithmetic_progression(self, lst: List[int]) -> bool:
        """
        Prüft ob eine Liste eine arithmetische Progression bildet.

        Eine AP hat konstante Differenz: lst[i+1] - lst[i] = konst.

        @param lst: Zu prüfende Liste
        @return: True falls lst eine AP ist
        @date: 2026-03-12
        """
        if len(lst) < 2:
            return True  # Trivial: 0 oder 1 Element ist immer AP
        d = lst[1] - lst[0]
        return all(lst[i + 1] - lst[i] == d for i in range(len(lst) - 1))

    def count_ap_free_colorings(self, N: int, k: int, r: int) -> int:
        """
        Zählt alle r-Färbungen von {1,...,N} ohne k-term AP (Brute-Force).

        ACHTUNG: Exponentiell in N*r — nur für kleine N praktikabel.

        @param N: Bereichsgröße (≤ 12 empfohlen)
        @param k: AP-Länge
        @param r: Anzahl Farben
        @return: Anzahl gültiger AP-freier Färbungen
        @date: 2026-03-12
        """
        count = 0
        # Alle r^N möglichen Färbungen durchprobieren
        for coloring in itertools.product(range(r), repeat=N):
            coloring_list = list(coloring)
            if self.is_ap_free_coloring(coloring_list, k):
                count += 1
        return count
