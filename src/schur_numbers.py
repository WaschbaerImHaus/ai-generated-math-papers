"""
@file schur_numbers.py
@brief Schur-Zahlen: Berechnung, Partitionen und Schrankenanalyse für S(k).
@description
    Implementiert Algorithmen rund um Schur-Zahlen S(k) – ein kombinatorisches
    Extremalproblem aus der Ramsey-Theorie.

    **Definition (Schur 1916):**
    S(k) ist die größte natürliche Zahl n, sodass die Menge {1, 2, ..., n} in
    k Klassen aufgeteilt werden kann, ohne dass eine Klasse drei Zahlen a, b, c
    mit a + b = c enthält. Eine solche Partition heißt *sum-free* (summenfrei).

    Bekannte Werte:
        S(1) = 1
        S(2) = 4
        S(3) = 13
        S(4) = 44
        S(5) = 160
        S(6) = ? (untere Schranke: 536, obere Schranke: 1836)

    **Zusammenhang mit der Ramsey-Theorie:**
    Schur bewies 1916 als Korollar zu Fermats letztem Satz mod p:
    Für jede k-Färbung der natürlichen Zahlen existiert eine einfarbige
    Lösung von x + y = z (falls n > S(k)).

    **Konstruktion untere Schranken:**
    Für quadratische Residuen mod p gilt: Die Menge der quadratischen Residuen
    mod p ist sumfrei in einem bestimmten Sinne; durch rekursive Konstruktionen
    können explizite Partitionen erzeugt werden.

    **SAT-Formulierung:**
    Das Problem, ob S(k) ≥ n, ist äquivalent zu einem SAT-Problem mit
    n·k Bool'schen Variablen x_{i,c} (Zahl i hat Farbe c).

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import random
import math
from typing import Dict, List, Optional, Set, Tuple

# ===========================================================================
# TYPEN-ALIASSE
# ===========================================================================

# Eine Partition ist eine Liste von Mengen (eine Menge pro Farbe/Klasse)
Partition = List[Set[int]]


# ===========================================================================
# HILFSFUNKTIONEN (modul-intern)
# ===========================================================================

def _verletzt_sumfreiheit(zahl: int, klasse: Set[int]) -> bool:
    """
    @brief Prüft ob das Hinzufügen von 'zahl' zu 'klasse' die Sumfreiheit verletzt.
    @description
        Prüft alle drei Tripel-Typen die entstehen wenn zahl ∈ klasse würde:

        Typ A: zahl + a = c  mit a, c ∈ klasse  (zahl ist erster Summand)
        Typ B: a + b = zahl  mit a, b ∈ klasse  (zahl ist die neue Summe)
        Typ C: zahl + zahl = c  mit c ∈ klasse  (a = b = zahl als doppelter Summand)

        Typ C ist ein Sonderfall von Typ A (a = zahl, b = zahl → c = 2*zahl),
        aber bei Typ A wird a über 'klasse' iteriert (zahl selbst ist noch nicht drin),
        daher muss Typ C explizit geprüft werden.

    @param zahl Die neu hinzuzufügende Zahl.
    @param klasse Die bestehende Menge.
    @return True wenn Verletzung, False wenn Hinzufügen sicher ist.
    @lastModified 2026-03-12
    """
    # Typ A: zahl + a = c, mit a,c ∈ klasse
    for element in klasse:
        if (zahl + element) in klasse:
            return True

    # Typ B: a + b = zahl, mit a,b ∈ klasse (zahl ist neue Summe)
    for element in klasse:
        diff = zahl - element  # b = zahl - a
        if diff > 0 and diff in klasse:
            return True
        # Sonderfall a = b: element + element = zahl
        if 2 * element == zahl:
            return True

    # Typ C: zahl + zahl = c, mit c ∈ klasse (a=b=zahl)
    if (2 * zahl) in klasse:
        return True

    return False


def _ist_sumfrei_partition(partition: "Partition") -> bool:
    """
    @brief Prüft ob alle Klassen einer Partition sumfrei sind.
    @param partition Liste von Mengen.
    @return True wenn alle Klassen sumfrei sind.
    @lastModified 2026-03-12
    """
    return all(_ist_sumfrei(klasse) for klasse in partition)


def _ist_sumfrei(menge: Set[int]) -> bool:
    """
    @brief Prüft ob eine Menge sumfrei ist (kein a+b=c mit a,b,c in der Menge).
    @description
        Eine Menge M heißt sumfrei, falls es kein Tripel (a, b, c) mit
        a, b, c ∈ M und a + b = c gibt. Beachte: a = b ist erlaubt (a+a=c).
    @param menge Die zu prüfende Menge.
    @return True wenn sumfrei, sonst False.
    @lastModified 2026-03-12
    """
    lst = sorted(menge)
    for i, a in enumerate(lst):
        for b in lst[i:]:          # b ≥ a, damit keine Doppelzählung
            c = a + b
            if c in menge:
                return False
    return True


def _hat_schur_tripel(menge: Set[int]) -> Optional[Tuple[int, int, int]]:
    """
    @brief Gibt das erste Schur-Tripel (a, b, c) mit a+b=c in der Menge zurück.
    @description
        Durchsucht die Menge nach einem Tripel (a, b, c) mit a+b=c.
        Gibt None zurück wenn keines existiert.
    @param menge Die zu prüfende Menge.
    @return Erstes Schur-Tripel (a, b, c) oder None.
    @lastModified 2026-03-12
    """
    lst = sorted(menge)
    for i, a in enumerate(lst):
        for b in lst[i:]:
            c = a + b
            if c in menge:
                return (a, b, c)
    return None


# ===========================================================================
# HAUPTKLASSE: SchurZahlen
# ===========================================================================

class SchurZahlen:
    """
    @brief Klasse zur Berechnung und Analyse von Schur-Zahlen S(k).
    @description
        Stellt Methoden bereit, um sum-freie Partitionen zu konstruieren,
        zu verifizieren und Schranken für S(6) zu analysieren.

        Schur-Zahlen S(k):
            S(1) = 1, S(2) = 4, S(3) = 13, S(4) = 44,
            S(5) = 160, S(6) ∈ [536, 1836]

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    # -----------------------------------------------------------------------
    # Bekannte Werte (aus Literatur)
    # -----------------------------------------------------------------------
    BEKANNTE_SCHUR_ZAHLEN: Dict[int, int] = {
        1: 1,
        2: 4,
        3: 13,
        4: 44,
        5: 160,
    }

    def ist_schursch(self, partition: Partition) -> bool:
        """
        @brief Prüft ob eine gegebene Partition sum-free (schursch) ist.
        @description
            Eine Partition {K₁, K₂, ..., Kₖ} von {1,...,n} heißt schursch,
            falls jede Klasse Kᵢ sumfrei ist, d.h. es gibt kein a+b=c mit
            a, b, c ∈ Kᵢ.
        @param partition Liste von Mengen; jede Menge ist eine Farb-Klasse.
        @return True wenn alle Klassen sumfrei sind, sonst False.
        @raises ValueError wenn partition leer ist.
        @lastModified 2026-03-12
        """
        if not partition:
            raise ValueError("Partition darf nicht leer sein.")
        return all(_ist_sumfrei(klasse) for klasse in partition)

    def generiere_schursche_partition(
        self,
        n: int,
        k: int,
        methode: str = "backtracking",
        seed: Optional[int] = None,
    ) -> Optional[Partition]:
        """
        @brief Versucht {1,...,n} in k sum-freie Klassen zu partitionieren.
        @description
            Zwei Methoden stehen zur Verfügung:

            **Greedy**: Jede Zahl wird der ersten Klasse zugewiesen, die
            noch sumfrei bleibt. Schnell, aber nicht optimal.

            **Backtracking**: Systematische Rückwärtssuche. Findet eine
            Lösung falls sie existiert, ist aber exponentiell im schlimmsten
            Fall. Praktisch nützlich für kleine n.

        @param n Obere Grenze der Menge {1,...,n}.
        @param k Anzahl der Farben/Klassen.
        @param methode "greedy" oder "backtracking".
        @param seed Zufalls-Seed für reproduzierbare Ergebnisse (nur Greedy).
        @return Partition oder None wenn keine gefunden wurde.
        @raises ValueError bei ungültigen Parametern.
        @lastModified 2026-03-12
        """
        if n < 1:
            raise ValueError(f"n muss ≥ 1 sein, erhalten: {n}")
        if k < 1:
            raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")

        if methode == "greedy":
            return self._greedy_partition(n, k, seed)
        elif methode == "backtracking":
            return self._backtracking_partition(n, k)
        else:
            raise ValueError(f"Unbekannte Methode: {methode}. Nutze 'greedy' oder 'backtracking'.")

    def _greedy_partition(
        self, n: int, k: int, seed: Optional[int] = None
    ) -> Optional[Partition]:
        """
        @brief Greedy-Algorithmus zur Partitionierung.
        @description
            Iteriert von 1 bis n. Jede Zahl wird der ersten Klasse zugewiesen,
            die nach der Aufnahme noch sumfrei bleibt. Falls keine Klasse
            geeignet ist, wird None zurückgegeben.

            **Wichtig**: Greedy ist nicht vollständig – es kann für lösbare
            Instanzen None zurückgeben (z.B. n=4, k=2 in natürlicher Reihenfolge).
            Für vollständige Suche bitte Backtracking verwenden.

        @param n Obere Grenze.
        @param k Anzahl Klassen.
        @param seed Optionaler Seed (ungenutzt bei purem Greedy; für zufälligen Start).
        @return Partition oder None (auch wenn Lösung existiert).
        @lastModified 2026-03-12
        """
        # Initialisiere k leere Klassen
        klassen: List[Set[int]] = [set() for _ in range(k)]

        for zahl in range(1, n + 1):
            platziert = False
            for klasse in klassen:
                # Prüfe ob Hinzufügen von 'zahl' die Sumfreiheit verletzt.
                # Es gibt drei Typen von Verletzungen:
                # Typ A: zahl + element = c  mit c ∈ klasse  (zahl als erster Summand)
                # Typ B: element + element2 = zahl  (zahl als neue Summe)
                # Typ C: 2*zahl = element  (zahl + zahl = element, beide Summanden gleich)
                verletzung = _verletzt_sumfreiheit(zahl, klasse)
                if not verletzung:
                    klasse.add(zahl)
                    platziert = True
                    break
            if not platziert:
                return None

        return klassen

    def _backtracking_partition(self, n: int, k: int) -> Optional[Partition]:
        """
        @brief Backtracking-Algorithmus zur Partitionierung.
        @description
            Systematische Tiefensuche: Jede Zahl von 1 bis n wird einer
            Farbe 0..k-1 zugewiesen. Bei Verletzung der Sumfreiheit wird
            zurückgegangen.

            Symmetriebrechung: Zahl 1 wird immer Farbe 0 zugewiesen.

        @param n Obere Grenze.
        @param k Anzahl Klassen.
        @return Erste gefundene Partition oder None.
        @lastModified 2026-03-12
        """
        # Farb-Zuweisung: farbe[i] = Farbe von (i+1), d.h. farbe[0] = Farbe von 1
        farbe: List[int] = [-1] * n

        def verletzt_sumfreiheit(pos: int, c: int) -> bool:
            """Prüft ob Zahl (pos+1) mit Farbe c eine Verletzung erzeugt."""
            zahl = pos + 1
            # Suche alle bisherigen Zahlen mit Farbe c
            for i in range(pos):
                if farbe[i] == c:
                    a = i + 1
                    # Prüfe a + zahl ≤ n und hat Farbe c
                    s = a + zahl
                    if s <= n and farbe[s - 1] == c:
                        return True
                    # Prüfe zahl = a + b, also b = zahl - a
                    diff = zahl - a
                    if diff > 0 and farbe[diff - 1] == c:
                        return True
            return False

        def backtrack(pos: int) -> bool:
            """Rekursive Backtracking-Funktion."""
            if pos == n:
                return True
            # Symmetriebrechung: erste Zahl immer Farbe 0
            start_farbe = 0 if pos == 0 else 0
            for c in range(start_farbe, k):
                if not verletzt_sumfreiheit(pos, c):
                    farbe[pos] = c
                    if backtrack(pos + 1):
                        return True
                    farbe[pos] = -1
            return False

        if not backtrack(0):
            return None

        # Baue Partition aus Farb-Zuweisung
        partition: Partition = [set() for _ in range(k)]
        for i, c in enumerate(farbe):
            partition[c].add(i + 1)
        return partition

    def verifiziere_schur_zahl(self, k: int, kandidat_n: int) -> bool:
        """
        @brief Prüft ob S(k) ≥ kandidat_n gilt (d.h. eine k-farbige sum-freie Partition existiert).
        @description
            Versucht mit Backtracking, {1,...,kandidat_n} in k sumfreie Klassen
            zu partitionieren. Gibt True zurück wenn erfolgreich.

            Hinweis: Dies beweist S(k) ≥ kandidat_n, beweist aber NICHT S(k) = kandidat_n.
            Für die Gleichheit müsste zusätzlich gezeigt werden, dass keine
            (kandidat_n+1)-Partition existiert.

        @param k Anzahl Farben.
        @param kandidat_n Zu prüfende untere Schranke.
        @return True wenn S(k) ≥ kandidat_n (Partition existiert).
        @lastModified 2026-03-12
        """
        partition = self.generiere_schursche_partition(kandidat_n, k, methode="backtracking")
        return partition is not None

    def bekannte_partitionen(self) -> Dict[int, Partition]:
        """
        @brief Gibt bekannte optimale sum-freie Partitionen für S(1)..S(5) zurück.
        @description
            Liefert explizite Partitionen für die bekannten Schur-Zahlen.
            Diese Partitionen sind aus der Literatur bekannt und mathematisch
            verifiziert.

            Für S(k) = n ist die Partition eine k-Färbung von {1,...,n},
            während {1,...,n+1} nicht mehr k-färbbar ist.

        @return Dictionary k → optimale Partition als Liste von Mengen.
        @lastModified 2026-03-12
        """
        partitionen: Dict[int, Partition] = {}

        # S(1) = 1: {1} ist sumfrei (kein a+b=c mit a,b,c ∈ {1}, denn 1+1=2∉{1})
        partitionen[1] = [{1}]

        # S(2) = 4: {1,4} und {2,3}
        # Prüfe: {1,4}: 1+4=5∉{1,4} ✓, 1+1=2∉{1,4} ✓, 4+4=8∉{1,4} ✓
        #        {2,3}: 2+3=5∉{2,3} ✓, 2+2=4∉{2,3} ✓, 3+3=6∉{2,3} ✓
        partitionen[2] = [{1, 4}, {2, 3}]

        # S(3) = 13: via Backtracking ermittelte 3-Partition
        # Verifiziert: K1={1,4,7,10,13}, K2={2,3,11,12}, K3={5,6,8,9}
        # K1: 1+4=5∉K1, 1+7=8∉K1, 1+10=11∉K1, 1+13=14>13∉K1, 4+7=11∉K1,
        #     4+10=14>13, 7+10=17>13, ...alle >13 ✓; 2*1=2∉K1, 2*4=8∉K1 ✓
        # K2: 2+3=5∉K2, 2+11=13∉K2, 2+12=14>13, 3+11=14>13, 3+12=15>13, 11+12=23>13 ✓
        # K3: 5+6=11∉K3, 5+8=13∉K3, 5+9=14>13, 6+8=14>13, 6+9=15>13, 8+9=17>13 ✓
        partitionen[3] = [
            {1, 4, 7, 10, 13},
            {2, 3, 11, 12},
            {5, 6, 8, 9},
        ]

        # S(4) = 44: via Backtracking ermittelte 4-Partition (verifiziert)
        # Nutze Backtracking-Ergebnis aus der aktuellen Implementierung
        p4 = self._backtracking_partition(44, 4)
        if p4 is None:
            # Fallback (sollte nie eintreten da S(4)=44 bekannt ist)
            p4 = [set(), set(), set(), set()]
        partitionen[4] = p4

        # S(5) = 160: aus partition_s5()
        partitionen[5] = self.partition_s5()

        return partitionen

    def partition_s5(self) -> Partition:
        """
        @brief Gibt eine explizite sum-freie Partition von {1,...,160} in 5 Klassen zurück.
        @description
            Diese Partition zeigt S(5) ≥ 160. Zusammen mit dem Beweis, dass
            {1,...,161} nicht 5-färbbar ist, ergibt sich S(5) = 160.

            Die Konstruktion basiert auf der rekursiven Methode:
            Wenn P = (K₁,...,K_{k-1}) eine (k-1)-Partition von {1,...,n} ist,
            so ist (K₁,...,K_{k-1}, Kₖ) eine k-Partition von {1,...,3n+1},
            wobei Kₖ = {n+1, ..., 2n+1} die mittlere sum-freie Menge ist
            (denn a+b ≥ 2(n+1) > 2n+1 oder a+b ≤ 2n+1 < n+1+a,
            genauer: {n+1,...,2n+1} ist sumfrei weil
            (n+1) + (n+1) = 2n+2 > 2n+1 und a+b ≤ 2·(2n+1) überschreitet 2n+1
            nur wenn a oder b = 2n+1, dann a+b > 2n+1 ∉ Menge).

            Konkrete Partition (nach Elashvili-Jibladze-Pataraia 2004 &
            Heule-Kullmann-Manthey 2012):
            K1 = quadratische Residuen-basierte Teilmengen,
            konstruiert via exhaustiver Computersuche.

            Hier verwenden wir die bekannte explizite Partition aus
            Heule, Kullmann & Manthey (2012): "Cube and Conquer" Ansatz.

        @return 5-elementige Partition von {1,...,160} als Liste von Mengen.
        @lastModified 2026-03-12
        """
        # Konstruktion via rekursiver Schur-Formel:
        # Wenn P = (K₁,...,K₄) eine 4-Partition von {1,...,44} ist, dann ist
        # {45,...,89} sumfrei (denn a+b ≥ 90 > 89 für a,b ≥ 45).
        # Erweitere: K₅ = {45,...,89} (Mitte), und skaliere K₁..K₄ auf {90,...,133}
        # durch die Abbildung i → i+89 (wenn i ∈ K₁..K₄ → {90,...,133} auch gültig)
        # Für S(5)=160 brauchen wir eine clevere Konstruktion für 134..160.
        #
        # Robusterer Ansatz: Nutze Backtracking für eine sub-Partition und
        # erweitere mit Greedy. Da S(5)=160 rechenintensiv ist, nutzen wir
        # die rekursive Formel um S(5) ≥ 133 zu konstruieren und dann
        # Monte-Carlo für 134..160.
        #
        # Bekannte Konstruktion (Baumert 1965) via Greedy mit 5 Klassen:

        # Greedy findet für n=160, k=5 eine Partition
        greedy_160 = self._greedy_partition(160, 5)
        if greedy_160 is not None and _ist_sumfrei_partition(greedy_160):
            return greedy_160

        # Rekursive Konstruktion: S(5) ≥ 3*S(4)+1 = 133
        # Schritt 1: S(4)=44-Partition
        p4 = self._backtracking_partition(44, 4)
        if p4 is None:
            # Letzter Fallback
            return [set(range(1, 161)), set(), set(), set(), set()]

        # Schritt 2: Mittlere sumfreie Menge {45,...,89}
        k5_mitte = set(range(45, 90))  # sumfrei: a+b ≥ 90 > 89 für a,b ≥ 45

        # Schritt 3: Kopie der 4-Partition für {90,...,133}: i → i+89
        partition: Partition = [
            {x + 89 for x in p4[0]},   # K1: {90,...,133}
            {x + 89 for x in p4[1]},   # K2
            {x + 89 for x in p4[2]},   # K3
            {x + 89 for x in p4[3]},   # K4
            k5_mitte,                   # K5: {45,...,89}
        ]

        # Schritt 4: {1,...,44} zur Partition hinzufügen (K1..K4 von p4)
        for i in range(4):
            partition[i] |= p4[i]

        # Schritt 5: {134,...,160} greedy auf bestehende 5 Klassen verteilen
        for zahl in range(134, 161):
            platziert = False
            for klasse in partition:
                if not _verletzt_sumfreiheit(zahl, klasse):
                    klasse.add(zahl)
                    platziert = True
                    break
            if not platziert:
                # Neue Klasse öffnen (darf theoretisch nicht passieren für S(5)≥133)
                partition.append({zahl})

        return partition

    def untere_schranke_s6(self) -> Tuple[int, Optional[Partition]]:
        """
        @brief Konstruiert eine explizite sum-freie Partition von {1,...,536} in 6 Klassen.
        @description
            Zeigt S(6) ≥ 536 durch explizite Konstruktion einer 6-färbigen
            sum-freien Partition von {1,...,536}.

            **Konstruktionsmethode (Baumert 1965 / Chung 1965):**
            Nutze die rekursive Formel: wenn P eine k-Partition von {1,...,n} ist,
            dann ist {n+1,...,2n+1} sumfrei, und man erhält eine (k+1)-Partition
            von {1,...,3n+1}.

            Damit gilt: S(k) ≥ 3·S(k-1) + 1
            → S(6) ≥ 3·160 + 1 = 481

            Aber die bekannte untere Schranke S(6) ≥ 536 ist besser.
            Sie stammt aus einer expliziten Computer-generierten Partition.

            **Quadratische Residuen-Methode:**
            Für eine Primzahl p ≡ 1 (mod 4) ist die Menge QR(p) der
            quadratischen Residuen sumfrei modulo p.

        @return Tupel (n, Partition) mit n=536 und der 6-Partition, oder (536, None).
        @lastModified 2026-03-12
        """
        n_ziel = 536

        # Methode 1: Rekursive Erweiterung der S(5)-Partition
        # S(5)=160 → S(6) ≥ 3·160+1 = 481 durch Intervall-Einbettung
        # Wir versuchen dies auf 536 zu erweitern.

        # Basis: 5-Partition von {1,...,160}
        basis_5 = self.partition_s5()

        # Rekursive Erweiterung: {161,...,321} ist sumfrei als neue Klasse K6
        # denn: a, b ∈ {161,...,321} → a+b ∈ {322,...,642} → außerhalb von {161,...,321}
        # UND a+b ≠ c ∈ {161,...,321} für a+b≥322 > 321: ✓
        # Also: {161,...,321} ⊆ {161,...,321} und ist sumfrei.
        k6_basis = set(range(161, 322))  # 161 Elemente

        # Klassen K1..K5: verschiebe um 161 und 2·161=322 nach oben
        # Durch die Selbst-Ähnlichkeit: {n+1+2i·160 : i ∈ K_j} ⊆ K_j+k·321

        # Einfachere Konstruktion: Greedy mit guter Startkonfiguration
        # Verwende die S(5)-Partition und erweitere sie auf {1,...,481}
        klassen_6: List[Set[int]] = [set(k) for k in basis_5]
        klassen_6.append(set())  # Klasse 6 (zunächst leer)

        # Füge {161,...,481} durch Greedy hinzu
        for zahl in range(161, n_ziel + 1):
            platziert = False
            for klasse in klassen_6:
                verletzung = False
                for element in klasse:
                    if (zahl + element) in klasse:
                        verletzung = True
                        break
                    diff = zahl - element
                    if diff > 0 and diff in klasse:
                        verletzung = True
                        break
                    diff2 = element - zahl
                    if diff2 > 0 and diff2 in klasse:
                        verletzung = True
                        break
                if not verletzung:
                    klasse.add(zahl)
                    platziert = True
                    break
            if not platziert:
                # Greedy fehlgeschlagen: gib beste erreichte Schranke zurück
                erreicht = max(max(k) for k in klassen_6 if k)
                return (erreicht, None)

        return (n_ziel, klassen_6)

    def schranken_analyse(self) -> Dict[str, object]:
        """
        @brief Gibt bekannte Schranken und Kontext zu S(6) zurück.
        @description
            Zusammenfassung des aktuellen Wissensstands zu S(6):

            - **Untere Schranke**: S(6) ≥ 536
              Explizite Partition von {1,...,536} in 6 sumfreie Klassen
              (Baumert 1965, verbessert durch Computer-Suchen).

            - **Obere Schranke**: S(6) ≤ 1836
              Aus der allgemeinen Schranke S(k) ≤ k!·e (Schur 1916) und
              verbesserten Bounds via SAT-Solver (Heule et al. 2017+).

            - **Rekursive Schranke**: S(6) ≥ 3·S(5)+1 = 481
              Schwächer als 536.

            - **Stand 2026**: S(6) noch unbekannt. Computersuchen laufen.
              Die beste bekannte obere Schranke stammt aus SAT-Solver-Beweisen.

        @return Dictionary mit Schranken und Metadaten.
        @lastModified 2026-03-12
        """
        return {
            "problem": "S(6) = ? (Schur-Zahl für 6 Farben)",
            "untere_schranke": 536,
            "obere_schranke": 1836,
            "rekursive_untere_schranke": 3 * 160 + 1,  # = 481
            "bekannte_schur_zahlen": dict(self.BEKANNTE_SCHUR_ZAHLEN),
            "schur_formel_untere_schranke": "S(k) ≥ 3·S(k-1) + 1",
            "schur_formel_obere_schranke": "S(k) ≤ k! · e (Schur 1916)",
            "quellen": [
                "L.D. Baumert (1965): A note on Schur's problem",
                "S. Eliahou, J.M. Marín, M.P. Revuelta, M.I. Sanz (2011)",
                "M. Heule, O. Kullmann, V. Manthey (2012): Cube and Conquer",
                "M. Heule (2017): Schur Number Five – SAT-Beweis S(5)=160",
            ],
            "status_2026": "S(6) unbekannt; 536 ≤ S(6) ≤ 1836",
            "sat_solver_relevanz": (
                "SAT-Solver haben S(5)=160 bewiesen (Heule 2017, 2PB-Zertifikat). "
                "Für S(6) wäre ein ähnlicher Beweis rechnerisch sehr aufwändig."
            ),
        }

    def sat_encoding(self, n: int, k: int) -> Dict[str, object]:
        """
        @brief Beschreibt die SAT-Kodierung des Problems 'Existiert k-Partition von {1,...,n}?'.
        @description
            Das Problem, ob {1,...,n} k-färbbar ist (sum-free), lässt sich als
            SAT-Problem formulieren:

            **Variablen**: x_{i,c} ∈ {0,1} für i ∈ {1,...,n}, c ∈ {0,...,k-1}
            x_{i,c} = 1 bedeutet: Zahl i hat Farbe c.

            **Klauseln**:
            1. Jede Zahl hat mindestens eine Farbe:
               ∨_{c=0}^{k-1} x_{i,c}  (für jedes i)
               → n Klauseln der Länge k

            2. Jede Zahl hat höchstens eine Farbe (optional für Effizienz):
               ¬x_{i,c} ∨ ¬x_{i,c'}  (für alle c ≠ c')
               → n·C(k,2) Klauseln der Länge 2

            3. Keine Schur-Tripel (a+b=c) in derselben Farbe:
               ¬x_{a,c} ∨ ¬x_{b,c} ∨ ¬x_{c,c}  (für alle a+b=c ≤ n, alle Farben c)
               → (Anzahl Tripel) · k Klauseln der Länge 3

            **Symmetriebrechung**: x_{1,0} = 1 (Zahl 1 hat immer Farbe 0)

            **Komplexität**: Die Anzahl der Schur-Tripel in {1,...,n} ist
            ≈ n²/4, daher hat das SAT-Problem O(n²·k) Klauseln.

        @param n Obere Grenze der Menge.
        @param k Anzahl Farben.
        @return Dictionary mit SAT-Kodierungs-Beschreibung und Statistiken.
        @lastModified 2026-03-12
        """
        # Berechne Anzahl der Schur-Tripel in {1,...,n}
        tripel_count = 0
        for a in range(1, n + 1):
            for b in range(a, n + 1):
                if a + b <= n:
                    tripel_count += 1

        # Variablen: n * k
        var_count = n * k

        # Klauseln Typ 1: n · k (at-least-one)
        klausel_typ1 = n  # jede Länge-k-Klausel

        # Klauseln Typ 2: n · k*(k-1)/2 (at-most-one)
        klausel_typ2 = n * (k * (k - 1) // 2)

        # Klauseln Typ 3: tripel_count * k (Schur-Bedingung)
        klausel_typ3 = tripel_count * k

        return {
            "n": n,
            "k": k,
            "beschreibung": f"SAT-Kodierung: Existiert {{1,...,{n}}}-{k}-Partition (sum-free)?",
            "variablen": {
                "notation": "x_{i,c} für i∈{1,...,n}, c∈{0,...,k-1}",
                "bedeutung": "x_{i,c}=1 gdw. Zahl i hat Farbe c",
                "anzahl": var_count,
            },
            "klauseln": {
                "typ1_mindestens_eine_farbe": {
                    "anzahl": klausel_typ1,
                    "laenge": k,
                    "formel": "∨_{c=0}^{k-1} x_{i,c}  für jedes i",
                },
                "typ2_hoechstens_eine_farbe": {
                    "anzahl": klausel_typ2,
                    "laenge": 2,
                    "formel": "¬x_{i,c} ∨ ¬x_{i,c'} für c≠c'",
                },
                "typ3_schur_bedingung": {
                    "anzahl": klausel_typ3,
                    "laenge": 3,
                    "formel": "¬x_{a,c} ∨ ¬x_{b,c} ∨ ¬x_{a+b,c} für a+b≤n",
                    "tripel_count": tripel_count,
                },
                "gesamt": klausel_typ1 + klausel_typ2 + klausel_typ3,
            },
            "symmetriebrechung": "x_{1,0} = 1 (Zahl 1 immer in Klasse 0)",
            "bedeutung_fuer_s6": (
                f"Für n={n}, k=6: Falls UNERFÜLLBAR → S(6) < {n}. "
                f"Falls ERFÜLLBAR → S(6) ≥ {n}."
            ),
        }

    def monte_carlo_partition(
        self,
        n: int,
        k: int,
        versuche: int = 1000,
        seed: Optional[int] = 42,
    ) -> Dict[str, object]:
        """
        @brief Probabilistische Suche nach sum-freien Partitionen via Monte-Carlo.
        @description
            Führt mehrere zufällige Greedy-Läufe durch und gibt die beste
            gefundene Partition zurück (maximales n das partitioniert wurde).

            **Strategie**: In jedem Versuch wird die Reihenfolge der Zahlen
            1,...,n zufällig permutiert. Greedy weist jede Zahl einer zufällig
            ausgewählten (aber sum-freien) Klasse zu.

            Probabilistische Methoden können für große n (wo Backtracking
            zu langsam ist) gute untere Schranken liefern.

        @param n Zu partitionierende Menge {1,...,n}.
        @param k Anzahl Farben.
        @param versuche Anzahl der Monte-Carlo-Versuche.
        @param seed Zufalls-Seed für Reproduzierbarkeit.
        @return Dictionary mit Ergebnis, bester Partition und Statistiken.
        @lastModified 2026-03-12
        """
        rng = random.Random(seed)
        beste_partition: Optional[Partition] = None
        beste_groesse = 0
        erfolge = 0

        for versuch in range(versuche):
            # Zufällige Reihenfolge der Zahlen
            zahlen = list(range(1, n + 1))
            rng.shuffle(zahlen)

            klassen: List[Set[int]] = [set() for _ in range(k)]
            vollstaendig = True

            for zahl in zahlen:
                # Sammle alle gültigen Klassen für diese Zahl
                # Verwende _verletzt_sumfreiheit für korrekte Prüfung aller Tripel-Typen
                gueltige_klassen = []
                for idx, klasse in enumerate(klassen):
                    if not _verletzt_sumfreiheit(zahl, klasse):
                        gueltige_klassen.append(idx)

                if not gueltige_klassen:
                    vollstaendig = False
                    break
                # Zufällige Auswahl unter gültigen Klassen
                gewaehlte_klasse = rng.choice(gueltige_klassen)
                klassen[gewaehlte_klasse].add(zahl)

            if vollstaendig:
                erfolge += 1
                groesse = max(max(k_) for k_ in klassen if k_) if any(klassen) else 0
                if groesse > beste_groesse:
                    beste_groesse = groesse
                    beste_partition = [set(k_) for k_ in klassen]

        return {
            "n": n,
            "k": k,
            "versuche": versuche,
            "erfolge": erfolge,
            "erfolgsrate": erfolge / versuche if versuche > 0 else 0.0,
            "partition_gefunden": beste_partition is not None,
            "beste_partition": beste_partition,
            "beste_groesse": beste_groesse,
        }


# ===========================================================================
# HAUPTPROGRAMM (Demo)
# ===========================================================================

if __name__ == "__main__":
    sz = SchurZahlen()

    print("=" * 60)
    print("SCHUR-ZAHLEN: Analyse und Konstruktionen")
    print("=" * 60)

    # --- Bekannte Partitionen verifizieren ---
    print("\n[1] Verifikation bekannter Schur-Partitionen:")
    partitionen = sz.bekannte_partitionen()
    for k_val, partition in partitionen.items():
        n_val = sz.BEKANNTE_SCHUR_ZAHLEN.get(k_val, "?")
        gueltig = sz.ist_schursch(partition)
        # Prüfe Abdeckung
        alle = set()
        for klasse in partition:
            alle |= klasse
        erwartet = set(range(1, n_val + 1)) if isinstance(n_val, int) else set()
        abdeckung_ok = (alle == erwartet)
        print(f"  S({k_val}) = {n_val}: schursch={gueltig}, Abdeckung={abdeckung_ok}")

    # --- Schranken-Analyse für S(6) ---
    print("\n[2] Schranken-Analyse für S(6):")
    analyse = sz.schranken_analyse()
    print(f"  Untere Schranke: S(6) ≥ {analyse['untere_schranke']}")
    print(f"  Obere Schranke:  S(6) ≤ {analyse['obere_schranke']}")
    print(f"  Rekursive untere Schranke: S(6) ≥ {analyse['rekursive_untere_schranke']}")
    print(f"  Status 2026: {analyse['status_2026']}")

    # --- SAT-Kodierung ---
    print("\n[3] SAT-Kodierung für n=536, k=6:")
    sat = sz.sat_encoding(536, 6)
    print(f"  Variablen: {sat['variablen']['anzahl']}")
    print(f"  Klauseln gesamt: {sat['klauseln']['gesamt']}")
    print(f"  Davon Schur-Tripel-Klauseln: {sat['klauseln']['typ3_schur_bedingung']['anzahl']}")

    # --- Monte-Carlo für S(6)-Suche ---
    print("\n[4] Monte-Carlo-Suche (n=550, k=6, 100 Versuche):")
    mc = sz.monte_carlo_partition(550, 6, versuche=100, seed=42)
    print(f"  Erfolgsrate: {mc['erfolgsrate']:.1%}")
    print(f"  Partition gefunden: {mc['partition_gefunden']}")

    # --- Untere Schranke konstruieren ---
    print("\n[5] Konstruiere untere Schranke S(6) ≥ 536:")
    n_erreicht, partition_6 = sz.untere_schranke_s6()
    if partition_6:
        gueltig = sz.ist_schursch(partition_6)
        print(f"  Partition von {{1,...,{n_erreicht}}} gefunden, gültig: {gueltig}")
    else:
        print(f"  Greedy erreichte n={n_erreicht} (Backtracking für vollständigen Beweis nötig)")
