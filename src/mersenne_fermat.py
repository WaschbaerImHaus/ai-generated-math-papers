"""
@file mersenne_fermat.py
@brief Mersenne-Primzahlen und Fermat-Primzahlen: Tests, Heuristiken und bekannte Ergebnisse.
@description
    Dieses Modul implementiert klassische Primzahltests und Strukturuntersuchungen
    für zwei besondere Primzahlklassen:

    **Mersenne-Primzahlen** (M_p = 2^p − 1):
    - Lucas-Lehmer-Test: Effizientester deterministischer Test für Mersenne-Primzahlen
    - Vollständige Liste aller 51 bekannten Mersenne-Primzahl-Exponenten
    - Wagstaff-Heuristik: Pr[M_p prim] ~ e^γ · ln(2) / ln(p)
    - Verbindung zu geraden vollkommenen Zahlen (Euklid-Euler-Theorem)

    **Fermat-Primzahlen** (F_n = 2^{2^n} + 1):
    - Pépin-Test: Deterministischer Primzahltest für Fermat-Zahlen
    - Nachweis: F_5 bis F_8 sind zusammengesetzt (bekannte Faktoren)
    - Heuristik: Wahrscheinlich nur F_0 bis F_4 prim (Conjecture)

    Mathematische Grundlagen:
    - Lucas-Lehmer: S_0 = 4, S_k = S_{k-1}^2 − 2 mod M_p; M_p prim ↔ S_{p-2} ≡ 0
    - Pépin: F_n prim ↔ 3^{(F_n−1)/2} ≡ −1 (mod F_n)
    - Euklid-Euler: n gerade vollkommen ↔ n = 2^{p-1}·(2^p − 1) mit M_p prim

@author Michael Fuhrmann
@lastModified 2026-03-12
"""

import math
from typing import Optional


# ──────────────────────────────────────────────────────────────────────────────
# Mersenne-Primzahlen
# ──────────────────────────────────────────────────────────────────────────────

class MersennePrimes:
    """
    @brief Klasse für Mersenne-Primzahlen M_p = 2^p − 1.

    @description
        Implementiert den Lucas-Lehmer-Test, stellt bekannte Exponenten bereit
        und berechnet Heuristiken sowie vollkommene Zahlen.

        **Mathematischer Hintergrund**:
        Eine Mersenne-Zahl M_p = 2^p − 1 kann nur dann prim sein, wenn p selbst
        prim ist (notwendige, aber nicht hinreichende Bedingung).

        Der Lucas-Lehmer-Test (1930/1878) ist der effizienteste bekannte Test
        für Mersenne-Zahlen und bildet die Grundlage aller GIMPS-Entdeckungen.

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    # Alle 51 bekannten Mersenne-Primzahl-Exponenten (Stand 2024)
    # Quelle: https://www.mersenne.org/primes/
    KNOWN_EXPONENTS = [
        2, 3, 5, 7, 13, 17, 19, 31, 61, 89,
        107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
        9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049,
        216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011,
        24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917,
        82589933
    ]

    # Eulersche Mascheroni-Konstante γ (benötigt für Wagstaff-Heuristik)
    EULER_MASCHERONI = 0.5772156649015328606

    def lucas_lehmer_test(self, p: int) -> bool:
        """
        @brief Führt den Lucas-Lehmer-Test für M_p = 2^p − 1 durch.

        @description
            **Algorithmus** (Lucas 1878, Lehmer 1930):
            1. Berechne M_p = 2^p − 1
            2. Setze S_0 = 4
            3. Iteriere: S_k = S_{k-1}^2 − 2 (mod M_p) für k = 1, ..., p−2
            4. M_p ist prim ↔ S_{p-2} ≡ 0 (mod M_p)

            **Satz** (Lucas-Lehmer, bewiesen):
            Für primes p ≥ 3 gilt: M_p = 2^p − 1 ist prim genau dann,
            wenn S_{p-2} ≡ 0 (mod M_p).

            Korrektheit: THEOREM (bewiesen). Kein Fehler möglich.

        @param p: Exponent (muss prim sein für nicht-triviale Aussage)
        @return True wenn M_p prim, sonst False
        @raises ValueError wenn p < 2
        """
        if p < 2:
            raise ValueError(f"Exponent p muss >= 2 sein, erhalten: {p}")

        # Sonderfall: M_2 = 3 ist prim (Test funktioniert auch, aber direkt klar)
        if p == 2:
            return True

        # p muss prim sein, sonst ist M_p zusammengesetzt
        if not self._is_prime(p):
            return False

        # Mersenne-Zahl berechnen
        mersenne = (1 << p) - 1  # 2^p − 1 (effizient mit Bitshift)

        # Lucas-Lehmer-Sequenz initialisieren
        s = 4

        # p−2 Iterationen (Indizes 0 bis p-3, also p-2 Schritte)
        for _ in range(p - 2):
            s = (s * s - 2) % mersenne

        return s == 0

    def mersenne_number(self, p: int) -> int:
        """
        @brief Berechnet die Mersenne-Zahl M_p = 2^p − 1.

        @param p: Exponent
        @return M_p als Integer
        """
        return (1 << p) - 1

    def is_known_mersenne_prime_exponent(self, p: int) -> bool:
        """
        @brief Prüft ob p ein bekannter Mersenne-Primzahl-Exponent ist.

        @param p: Zu prüfender Exponent
        @return True wenn p in der Liste der 51 bekannten Exponenten
        """
        return p in self.KNOWN_EXPONENTS

    def wagstaff_heuristic(self, p: int) -> float:
        """
        @brief Berechnet die Wagstaff-Heuristik für die Primzahlwahrscheinlichkeit von M_p.

        @description
            **Wagstaff-Heuristik** (Conjecture, nicht bewiesen):
            Die Wahrscheinlichkeit, dass M_p = 2^p − 1 eine Primzahl ist,
            beträgt näherungsweise:

            Pr[M_p prim] ≈ e^γ · ln(2) / ln(p)

            wobei γ ≈ 0.5772... die Euler-Mascheroni-Konstante ist.

            **ACHTUNG**: Dies ist eine CONJECTURE (Wagstaff 1983), kein bewiesener Satz.
            Die Heuristik basiert auf probabilistischen Modellen der Primzahlverteilung
            und dem Satz von Mertens, ist aber nicht rigoros bewiesen.

            Implikation: Unendlich viele Mersenne-Primzahlen (Conjecture!).

        @param p: Primexponent
        @return Geschätzte Wahrscheinlichkeit (zwischen 0 und 1)
        @raises ValueError wenn p <= 1
        """
        if p <= 1:
            raise ValueError(f"p muss > 1 sein, erhalten: {p}")

        # e^γ · ln(2) / ln(p)
        e_gamma = math.exp(self.EULER_MASCHERONI)
        probability = e_gamma * math.log(2) / math.log(p)

        # Wahrscheinlichkeit ist maximal 1
        return min(probability, 1.0)

    def perfect_number_from_mersenne(self, p: int) -> Optional[int]:
        """
        @brief Berechnet die gerade vollkommene Zahl zu M_p (falls M_p prim).

        @description
            **Euklid-Euler-Theorem** (THEOREM, vollständig bewiesen):
            Eine gerade Zahl n ist vollkommen genau dann, wenn gilt:
            n = 2^{p-1} · (2^p − 1) mit M_p = 2^p − 1 prim.

            - Euklid (~300 v. Chr.): Wenn M_p prim, dann ist 2^{p-1}·M_p vollkommen.
            - Euler (1747): Jede gerade vollkommene Zahl hat diese Form.

            Beispiele:
            - p=2: 2^1 · 3 = 6 (kleinste vollkommene Zahl)
            - p=3: 2^2 · 7 = 28
            - p=5: 2^4 · 31 = 496
            - p=7: 2^6 · 127 = 8128

            Ob ungerade vollkommene Zahlen existieren, ist UNBEKANNT (offenes Problem).

        @param p: Exponent (M_p wird auf Primalität getestet)
        @return Vollkommene Zahl wenn M_p prim, sonst None
        """
        if self.lucas_lehmer_test(p):
            mersenne = self.mersenne_number(p)
            return (1 << (p - 1)) * mersenne  # 2^{p-1} · M_p
        return None

    def verify_perfect_number(self, n: int) -> bool:
        """
        @brief Verifiziert ob n eine vollkommene Zahl ist (σ(n) = 2n).

        @description
            Eine Zahl n heißt vollkommen, wenn die Summe aller echten Teiler
            gleich n ist, äquivalent: σ(n) = 2n.

        @param n: Zu prüfende Zahl
        @return True wenn n vollkommen
        """
        if n <= 0:
            return False

        # Summe aller Teiler einschließlich 1, aber ohne n selbst
        divisor_sum = sum(d for d in range(1, n) if n % d == 0)
        return divisor_sum == n

    def get_first_n_mersenne_primes(self, n: int) -> list:
        """
        @brief Gibt die ersten n Mersenne-Primzahl-Exponenten zurück.

        @param n: Anzahl der gewünschten Exponenten
        @return Liste der ersten n bekannten Mersenne-Primzahl-Exponenten
        @raises ValueError wenn n > 51 (nur 51 bekannte)
        """
        if n > len(self.KNOWN_EXPONENTS):
            raise ValueError(
                f"Nur {len(self.KNOWN_EXPONENTS)} Mersenne-Primzahlen bekannt, "
                f"aber {n} angefordert"
            )
        return self.KNOWN_EXPONENTS[:n]

    def digit_count(self, p: int) -> int:
        """
        @brief Berechnet die Anzahl der Dezimalstellen von M_p = 2^p − 1.

        @description
            Anzahl Stellen = floor(p · log10(2)) + 1

        @param p: Exponent
        @return Anzahl Dezimalstellen
        """
        return int(p * math.log10(2)) + 1

    def _is_prime(self, n: int) -> bool:
        """
        @brief Einfacher Primzahltest (Trial Division) für interne Verwendung.

        @param n: Zu prüfende Zahl
        @return True wenn n prim
        """
        if n < 2:
            return False
        if n == 2:
            return True
        if n % 2 == 0:
            return False
        for i in range(3, int(math.isqrt(n)) + 1, 2):
            if n % i == 0:
                return False
        return True


# ──────────────────────────────────────────────────────────────────────────────
# Fermat-Primzahlen
# ──────────────────────────────────────────────────────────────────────────────

class FermatPrimes:
    """
    @brief Klasse für Fermat-Primzahlen F_n = 2^{2^n} + 1.

    @description
        Fermat-Zahlen haben die Form F_n = 2^{2^n} + 1:
        - F_0 = 3 (prim)
        - F_1 = 5 (prim)
        - F_2 = 17 (prim)
        - F_3 = 257 (prim)
        - F_4 = 65537 (prim)
        - F_5 = 4294967297 = 641 · 6700417 (zusammengesetzt, Euler 1732)
        - F_6–F_32: alle bekannt zusammengesetzt

        **Fermat-Conjecture** (UNBEWIESEN):
        Es gibt nur endlich viele Fermat-Primzahlen. Wahrscheinlich sind
        F_0 bis F_4 die einzigen (Heuristik: Σ Pr[F_n prim] konvergiert).

        Verbindung zur Geometrie: F_n prim ↔ reguläres 2^n·F_n-Eck konstruierbar
        (Gauß-Wantzel-Theorem, BEWIESEN).

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    # Bekannte Faktoren zusammengesetzter Fermat-Zahlen
    # Format: {n: [(faktor1, exponent1), ...]}
    KNOWN_FACTORS = {
        5:  [(641, 1), (6700417, 1)],           # Euler 1732
        6:  [(274177, 1), (67280421310721, 1)],  # Landry & Le Lasseur 1880
        7:  [(59649589127497217, 1), (5704689200685129054721, 1)],  # Morrison & Brillhart 1975
        8:  [(1238926361552897, 1)],              # Brent & Pollard 1981 (teilfaktorisiert)
    }

    # Bekannte Primzahl-Fermat-Zahlen (VOLLSTÄNDIG BEWIESEN für n=0..4)
    KNOWN_PRIME_INDICES = [0, 1, 2, 3, 4]

    def fermat_number(self, n: int) -> int:
        """
        @brief Berechnet die n-te Fermat-Zahl F_n = 2^{2^n} + 1.

        @param n: Index (n >= 0)
        @return F_n als Integer
        @raises ValueError wenn n < 0
        """
        if n < 0:
            raise ValueError(f"Fermat-Index n muss >= 0 sein, erhalten: {n}")

        # 2^{2^n} + 1 (für n >= 20 wird das sehr groß!)
        return (1 << (1 << n)) + 1

    def pepin_test(self, n: int) -> bool:
        """
        @brief Führt den Pépin-Test für F_n = 2^{2^n} + 1 durch.

        @description
            **Pépin-Test** (THEOREM, Pépin 1877):
            Für n ≥ 1 gilt: F_n = 2^{2^n} + 1 ist prim genau dann, wenn

            3^{(F_n − 1)/2} ≡ −1 (mod F_n)

            Beweis basiert auf dem Euler-Kriterium: Für eine Primzahl p und
            a mit gcd(a,p)=1 gilt: a^{(p-1)/2} ≡ (a/p) (mod p),
            wobei (a/p) das Legendre-Symbol ist.

            Da 3 ein quadratisches Nichtrest modulo F_n ist (für n≥1, sofern F_n prim),
            ist das Legendre-Symbol (3/F_n) = −1, was den Test ergibt.

            **Effizienz**: (F_n − 1)/2 = 2^{2^n − 1}, also sind 2^n − 1 modulare
            Quadrierungen nötig (schnelle Exponentiation mit pow()).

            **ACHTUNG**: Für n ≥ 5 sind die Zahlen astronomisch groß (F_5 hat 10 Stellen,
            F_20 hat über 300.000 Stellen). Daher max_n = 8 als Grenze.

        @param n: Fermat-Index
        @return True wenn F_n prim, sonst False
        @raises ValueError wenn n < 0 oder n > 8 (zu groß für praktische Berechnung)
        """
        if n < 0:
            raise ValueError(f"n muss >= 0 sein, erhalten: {n}")
        if n > 8:
            raise ValueError(
                f"n={n} zu groß für Pépin-Test (F_{n} hat astronomisch viele Stellen). "
                f"Bekannte Faktoren verwenden."
            )

        # Sonderfall F_0 = 3 (prim, Pépin-Test gilt strikt für n >= 1)
        if n == 0:
            return True  # F_0 = 3 ist offensichtlich prim

        fn = self.fermat_number(n)

        # Exponent: (F_n − 1) / 2 = 2^{2^n − 1}
        exponent = fn >> 1  # (F_n − 1) / 2, denn F_n − 1 = 2^{2^n} ist gerade

        # Modulare Exponentiation: 3^exponent mod F_n
        result = pow(3, exponent, fn)

        # Prim ↔ Ergebnis ≡ −1 ≡ F_n − 1 (mod F_n)
        return result == fn - 1

    def is_prime_fermat(self, n: int) -> bool:
        """
        @brief Prüft ob F_n prim ist (kombiniert bekannte Ergebnisse + Pépin-Test).

        @description
            Verwendet zunächst bekannte Faktoren, dann den Pépin-Test.

        @param n: Fermat-Index
        @return True wenn F_n prim
        """
        # Bekannte Primzahl-Indices
        if n in self.KNOWN_PRIME_INDICES:
            return True

        # Bekannte zusammengesetzte Zahlen
        if n in self.KNOWN_FACTORS:
            return False

        # Pépin-Test für kleinere n
        if n <= 8:
            return self.pepin_test(n)

        # Für n > 8: keine bekannte Entscheidung ohne Spezialhardware
        raise NotImplementedError(
            f"F_{n} ist zu groß für diesen Test. Verwende spezialisierte Software."
        )

    def get_known_factor(self, n: int) -> Optional[int]:
        """
        @brief Gibt einen bekannten Faktor von F_n zurück (falls vorhanden).

        @param n: Fermat-Index
        @return Kleinster bekannter Faktor oder None wenn unbekannt/prim
        """
        if n in self.KNOWN_FACTORS:
            return self.KNOWN_FACTORS[n][0][0]  # Erster (kleinster) Faktor
        return None

    def verify_factorization(self, n: int) -> bool:
        """
        @brief Verifiziert die bekannte Faktorisierung von F_n.

        @description
            Prüft, ob der bekannte Faktor tatsächlich F_n teilt.

        @param n: Fermat-Index (muss in KNOWN_FACTORS sein)
        @return True wenn Faktor korrekt
        @raises KeyError wenn n nicht in KNOWN_FACTORS
        """
        if n not in self.KNOWN_FACTORS:
            raise KeyError(f"Keine bekannten Faktoren für F_{n}")

        fn = self.fermat_number(n)
        factor = self.KNOWN_FACTORS[n][0][0]

        return fn % factor == 0

    def wagstaff_heuristic_fermat(self) -> dict:
        """
        @brief Berechnet Wagstaff-Heuristik: Warum wahrscheinlich nur F_0−F_4 prim.

        @description
            **Heuristik** (CONJECTURE, nicht bewiesen):
            Die Wahrscheinlichkeit, dass F_n prim ist, beträgt näherungsweise:

            Pr[F_n prim] ≈ 1 / (2^n · ln(2))

            (aus dem Primzahlsatz: π(x) ~ x/ln(x), angewendet auf Fermat-Zahlen)

            Die Summe Σ_{n=0}^∞ Pr[F_n prim] konvergiert, was suggeriert,
            dass nur endlich viele Fermat-Primzahlen existieren.

            **WICHTIG**: Dies ist eine HEURISTIK, kein Beweis. Es ist nicht
            bewiesen, dass F_5, F_6, ... alle zusammengesetzt sind (obwohl
            alle bisher untersuchten Fermat-Zahlen F_5 bis F_32 zusammengesetzt sind).

        @return Dictionary mit heuristischen Wahrscheinlichkeiten für n=0..20
        """
        result = {}
        for n in range(21):
            # Pr[F_n prim] ≈ 1 / (2^n · ln(2))
            if n == 0:
                prob = 1.0  # F_0 = 3 ist definitiv prim
            else:
                prob = 1.0 / ((1 << n) * math.log(2))
            result[n] = min(prob, 1.0)

        # Gesamtwahrscheinlichkeit (Summe = endlich → endlich viele Primzahlen zu erwarten)
        result['total_expected'] = sum(v for k, v in result.items() if isinstance(k, int))
        return result

    def gauss_wantzel_constructible(self, n: int) -> bool:
        """
        @brief Prüft ob ein reguläres 2^n-Eck mit Zirkel+Lineal konstruierbar ist.

        @description
            **Gauß-Wantzel-Theorem** (BEWIESEN, Gauß 1796 / Wantzel 1837):
            Ein reguläres p-Eck ist mit Zirkel und Lineal konstruierbar genau dann,
            wenn p ein Produkt aus einer Zweierpotenz und verschiedenen Fermat-Primzahlen ist.

            Insbesondere: Reguläres F_n-Eck konstruierbar ↔ F_n ist prim.

        @param n: Fermat-Index
        @return True wenn F_n-Eck konstruierbar (d.h. wenn F_n prim ist)
        """
        if n in self.KNOWN_PRIME_INDICES:
            return True
        if n in self.KNOWN_FACTORS:
            return False
        # Für unbekannte F_n: Conjecture = False (alle F_n mit n>=5 vermutlich zusammengesetzt)
        return False  # Conjecture, nicht bewiesen
