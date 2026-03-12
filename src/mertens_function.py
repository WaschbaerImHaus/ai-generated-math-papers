"""
Mertens-Funktion M(x), Liouville-Funktion λ(n) und verwandte Funktionen.

Dieses Modul implementiert:
    - Möbius-Funktion μ(n) via Siebmethode
    - Mertens-Funktion M(x) = Σ_{n≤x} μ(n)
    - Liouville-Funktion λ(n) = (−1)^Ω(n)
    - Summatorische Liouville-Funktion L(x) = Σ_{n≤x} λ(n)
    - Pólya-Vermutung (widerlegt 1960 von Haselgrove)
    - Mertens-Vermutung |M(x)| < √x (widerlegt 1985 von Odlyzko/te Riele)
    - Schranken unter Riemann-Hypothese: M(x) = O(x^{1/2+ε})

Mathematischer Hintergrund:
    Die Mertens-Funktion ist eng mit der Riemann-Hypothese verbunden:
    RH ⟺ M(x) = O(x^{1/2+ε}) für alle ε > 0.
    Die stärkere Mertens-Vermutung |M(x)| < √x wurde 1985 von Odlyzko und
    te Riele widerlegt – allerdings existiert bis heute kein explizites Gegenbeispiel.
    Das kleinste ist (heuristisch) bei x ≈ exp(1.59·10^{40}).

    Pintz (1987) zeigte: lim-sup_{x→∞} M(x)/√x > 1.06...

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import math
import numpy as np
from typing import List, Optional, Tuple, Dict
from functools import lru_cache


class MertensFunction:
    """
    Berechnet die Mertens-Funktion M(x) = Σ_{n≤x} μ(n).

    Die Möbius-Funktion μ(n) ist definiert als:
        μ(1) = 1
        μ(n) = (−1)^k falls n = p₁·p₂·...·pₖ (quadratfrei mit k Primfaktoren)
        μ(n) = 0 falls n durch ein Primquadrat teilbar ist

    Die Mertens-Funktion akkumuliert μ:
        M(x) = Σ_{n=1}^{⌊x⌋} μ(n)

    Bekannte Werte:
        M(10)  = −1
        M(100) = 1
        M(1000) = 2

    CONJECTURE (Mertens 1897, widerlegt 1985):
        |M(x)| < √x für alle x ≥ 1 → FALSCH (Odlyzko/te Riele)

    THEOREM (unter RH):
        M(x) = O(x^{1/2+ε}) für alle ε > 0

    THEOREM (Pintz 1987):
        lim-sup_{x→∞} M(x)/√x > 1.06

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self, limit: int):
        """
        Initialisiert den Sieb bis zum Limit und berechnet μ(n) für alle n ≤ limit.

        Der lineare Sieb (Sieb des Eratosthenes für μ) hat Zeitkomplexität O(N log log N)
        und Speicherkomplexität O(N).

        @param limit: Obere Grenze für die Siebberechnung (inklusive)
        @lastModified: 2026-03-12
        """
        if limit < 1:
            raise ValueError("Limit muss mindestens 1 sein.")
        self.limit = limit
        # Sieb: Berechne μ(n) für n = 1..limit
        self._mu: np.ndarray = self._sieve_mobius(limit)
        # Kumulative Summe: M(x) = Σ_{n=1}^{x} μ(n)
        # self._M[i] enthält M(i) für i = 0,1,...,limit
        self._M: np.ndarray = np.cumsum(self._mu)

    @staticmethod
    def _sieve_mobius(n: int) -> np.ndarray:
        """
        Berechnet μ(k) für k = 1, ..., n via linearem Sieb.

        Algorithmus:
            1. Initialisiere smallest_prime[k] = 0
            2. Durchlaufe Primzahlen p via Sieb
            3. Für jede Zahl k: teile durch kleinsten Primfaktor, prüfe Potenz

        Rückgabe: Array mu mit mu[k-1] = μ(k) für k = 1..n

        @param n: Obergrenze
        @return: numpy-Array der Länge n, Index 0 entspricht μ(1)
        @lastModified: 2026-03-12
        """
        # mu[i] = μ(i+1), also mu[0] = μ(1) = 1
        mu = np.zeros(n + 1, dtype=np.int8)
        mu[1] = 1

        # Zähle kleinste Primfaktoren-Potenzen für jeden Wert
        # is_composite[i] = True falls i zusammengesetzt
        is_square_free = np.ones(n + 1, dtype=bool)  # Annahme: quadratfrei
        primes = []
        spf = np.arange(n + 1, dtype=np.int32)  # kleinster Primfaktor (smallest prime factor)

        # Linearer Sieb
        for i in range(2, n + 1):
            if spf[i] == i:
                # i ist Primzahl
                primes.append(i)
                mu[i] = -1  # μ(p) = -1 für Primzahl p

            for p in primes:
                if p > spf[i] or i * p > n:
                    break
                spf[i * p] = p
                if i % p == 0:
                    # p^2 | i*p → nicht quadratfrei
                    mu[i * p] = 0
                else:
                    # μ ist multiplikativ: μ(i*p) = μ(i)*μ(p) = μ(i)*(-1)
                    mu[i * p] = -mu[i]

        # Rückgabe: mu[1..n], Index 0 auf 0 (ungenutzt)
        return mu[1:]  # mu[0] = μ(1), ..., mu[n-1] = μ(n)

    def mobius(self, n: int) -> int:
        """
        Gibt μ(n) zurück.

        @param n: Positive ganze Zahl ≤ self.limit
        @return: Möbius-Wert ∈ {-1, 0, 1}
        @lastModified: 2026-03-12
        """
        if n < 1 or n > self.limit:
            raise ValueError(f"n={n} außerhalb des Bereichs [1, {self.limit}]")
        return int(self._mu[n - 1])

    def M(self, x: float) -> int:
        """
        Berechnet M(x) = Σ_{n=1}^{⌊x⌋} μ(n).

        @param x: Reelle Zahl ≥ 1
        @return: Ganzzahliger Wert der Mertens-Funktion
        @lastModified: 2026-03-12
        """
        k = int(math.floor(x))
        if k < 1:
            return 0
        if k > self.limit:
            raise ValueError(f"x={x} > limit={self.limit}. Bitte größeres Limit wählen.")
        return int(self._M[k - 1])

    def M_array(self, x_max: Optional[int] = None) -> np.ndarray:
        """
        Gibt das Array [M(1), M(2), ..., M(x_max)] zurück.

        Nützlich für Plotfunktionen und statistische Analysen.

        @param x_max: Maximaler x-Wert (Standard: self.limit)
        @return: Array der Mertens-Werte
        @lastModified: 2026-03-12
        """
        if x_max is None:
            x_max = self.limit
        return self._M[:x_max].copy()

    def ratio_M_sqrt(self, x_max: Optional[int] = None) -> np.ndarray:
        """
        Berechnet M(x)/√x für x = 1, ..., x_max.

        Dient der visuellen Prüfung der (widerlegten) Mertens-Vermutung:
        CONJECTURE (FALSCH): |M(x)/√x| < 1 für alle x

        @param x_max: Maximaler x-Wert
        @return: Array von M(x)/√x
        @lastModified: 2026-03-12
        """
        if x_max is None:
            x_max = self.limit
        x_vals = np.arange(1, x_max + 1, dtype=float)
        m_vals = self._M[:x_max].astype(float)
        return m_vals / np.sqrt(x_vals)

    def max_ratio(self, x_max: Optional[int] = None) -> Tuple[int, float]:
        """
        Findet das x ≤ x_max, wo |M(x)/√x| maximal ist.

        Relevant für Odlyzko/te Riele Widerlegung der Mertens-Vermutung.

        @param x_max: Maximaler x-Wert
        @return: Tupel (x_max_ratio, max_ratio_value)
        @lastModified: 2026-03-12
        """
        if x_max is None:
            x_max = self.limit
        ratios = np.abs(self.ratio_M_sqrt(x_max))
        idx = int(np.argmax(ratios))
        return idx + 1, float(ratios[idx])

    def verify_mertens_conjecture(self, x_max: Optional[int] = None) -> bool:
        """
        Prüft, ob |M(x)| < √x für alle x ≤ x_max gilt.

        HINWEIS: Die Mertens-Vermutung ist für x < ≈ 10^{13} numerisch gültig,
        aber theoretisch widerlegt (Odlyzko/te Riele 1985).

        @param x_max: Obere Grenze der Überprüfung
        @return: True wenn die Ungleichung für alle geprüften x gilt
        @lastModified: 2026-03-12
        """
        if x_max is None:
            x_max = self.limit
        x_vals = np.arange(1, x_max + 1, dtype=float)
        m_vals = np.abs(self._M[:x_max].astype(float))
        return bool(np.all(m_vals < np.sqrt(x_vals)))

    def zeros_and_M_connection(self) -> str:
        """
        Erklärt die theoretische Verbindung zwischen Mertens-Funktion und RH-Nullstellen.

        Explizite Formel (formal, unter RH):
            M(x) ~ Σ_{ρ: ζ(ρ)=0} x^ρ / (ρ · ζ'(ρ))

        Diese Formel zeigt: Wachstum von M(x) ist durch Real-Teil der Nullstellen kontrolliert.
        Unter RH: Re(ρ) = 1/2 → M(x) = O(x^{1/2+ε})

        @return: Erklärungstext
        @lastModified: 2026-03-12
        """
        return (
            "Explicit formula: M(x) = Σ_{ρ} x^ρ/(ρ·ζ'(ρ)) − 2 + ...\n"
            "Unter RH (Re(ρ)=1/2): M(x) = O(x^{1/2+ε}) für alle ε>0\n"
            "Pintz (1987): lim-sup M(x)/√x > 1.06 (THEOREM)\n"
            "Odlyzko/te Riele (1985): Mertens-Vermutung WIDERLEGT\n"
            "Kein explizites Gegenbeispiel bekannt (heuristisch: x ≈ exp(1.59·10^40))"
        )


class LiouvilleFunction:
    """
    Liouville-Funktion λ(n) = (−1)^Ω(n) und summatorische Funktion L(x).

    Ω(n) = Gesamtzahl der Primfaktoren von n (mit Vielfachheit).

    Beispiele:
        λ(1) = 1    (Ω(1) = 0, leeres Produkt)
        λ(2) = -1   (Ω(2) = 1, 2 = 2)
        λ(4) = 1    (Ω(4) = 2, 4 = 2²)
        λ(6) = 1    (Ω(6) = 2, 6 = 2·3)
        λ(12) = 1   (Ω(12) = 3, 12 = 2²·3)

    Summatorische Funktion: L(x) = Σ_{n=1}^{⌊x⌋} λ(n)

    CONJECTURE (Pólya 1919, widerlegt 1960):
        L(x) ≤ 0 für alle x ≥ 2 → FALSCH (Haselgrove 1960)
        Kleinstes bekanntes Gegenbeispiel: x ≈ 906316571

    THEOREM (Liouville, Tschebyschow):
        Σ_{n≤x} λ(n)/n = Ο(1) (äquivalent zu PNT)

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self, limit: int):
        """
        Initialisiert den Sieb für Liouville-Funktion bis zum Limit.

        @param limit: Obergrenze (inklusive)
        @lastModified: 2026-03-12
        """
        if limit < 1:
            raise ValueError("Limit muss mindestens 1 sein.")
        self.limit = limit
        self._lambda: np.ndarray = self._sieve_liouville(limit)
        self._L: np.ndarray = np.cumsum(self._lambda)

    @staticmethod
    def _sieve_liouville(n: int) -> np.ndarray:
        """
        Berechnet λ(k) = (−1)^Ω(k) für k = 1..n via Sieb.

        Algorithmus: Für jede Zahl k, dividiere wiederholt durch kleinsten
        Primfaktor und zähle Primfaktoren (mit Vielfachheit).

        @param n: Obergrenze
        @return: Array der Länge n, Index i entspricht λ(i+1)
        @lastModified: 2026-03-12
        """
        # omega[k] = Ω(k) = Anzahl Primfaktoren mit Vielfachheit
        omega = np.zeros(n + 1, dtype=np.int32)

        # Sieb für Ω(n): Addiere 1 für jeden Primfaktor-Treffer
        for p in range(2, n + 1):
            if omega[p] == 0:
                # p ist Primzahl (noch nicht bearbeitet)
                # Erhöhe Ω für alle Vielfachen von p^k
                pk = p
                while pk <= n:
                    # Für alle Vielfachen von p^k: +1 zum Ω-Zähler
                    for mult in range(pk, n + 1, pk):
                        omega[mult] += 1
                    pk *= p  # nächste Potenz

        # λ(k) = (-1)^Ω(k)
        lam = np.where(omega % 2 == 0, 1, -1).astype(np.int8)
        lam[0] = 0  # Index 0 ungenutzt
        return lam[1:]  # lam[0] = λ(1), ..., lam[n-1] = λ(n)

    def liouville(self, n: int) -> int:
        """
        Gibt λ(n) zurück.

        @param n: Positive ganze Zahl ≤ self.limit
        @return: λ(n) ∈ {-1, 1}
        @lastModified: 2026-03-12
        """
        if n < 1 or n > self.limit:
            raise ValueError(f"n={n} außerhalb [1, {self.limit}]")
        return int(self._lambda[n - 1])

    def big_omega(self, n: int) -> int:
        """
        Berechnet Ω(n) = Gesamtzahl der Primfaktoren (mit Vielfachheit).

        Beispiel: Ω(12) = 3, da 12 = 2·2·3.

        @param n: Positive ganze Zahl
        @return: Ω(n)
        @lastModified: 2026-03-12
        """
        if n < 1:
            raise ValueError("n muss positiv sein")
        count = 0
        k = n
        d = 2
        while d * d <= k:
            while k % d == 0:
                count += 1
                k //= d
            d += 1
        if k > 1:
            count += 1
        return count

    def L(self, x: float) -> int:
        """
        Berechnet L(x) = Σ_{n=1}^{⌊x⌋} λ(n).

        @param x: Reelle Zahl ≥ 1
        @return: Summatorische Liouville-Funktion
        @lastModified: 2026-03-12
        """
        k = int(math.floor(x))
        if k < 1:
            return 0
        if k > self.limit:
            raise ValueError(f"x > limit={self.limit}")
        return int(self._L[k - 1])

    def L_array(self, x_max: Optional[int] = None) -> np.ndarray:
        """
        Gibt [L(1), L(2), ..., L(x_max)] zurück.

        @param x_max: Maximaler x-Wert
        @return: Array der L-Werte
        @lastModified: 2026-03-12
        """
        if x_max is None:
            x_max = self.limit
        return self._L[:x_max].copy()

    def polya_conjecture_check(self, x_max: Optional[int] = None) -> Dict[str, object]:
        """
        Prüft die (widerlegte) Pólya-Vermutung L(x) ≤ 0 für x ≥ 2.

        CONJECTURE (Pólya 1919, WIDERLEGT 1960 von Haselgrove):
            L(x) ≤ 0 für alle x ≥ 2

        Kleinstes bekanntes Gegenbeispiel: x ≈ 906316571

        @param x_max: Obere Grenze der Überprüfung
        @return: Dict mit Status und ggf. Gegenbeispiel
        @lastModified: 2026-03-12
        """
        if x_max is None:
            x_max = self.limit

        # Suche nach positiven L(x) für x >= 2
        violations = []
        for x in range(2, x_max + 1):
            val = int(self._L[x - 1])
            if val > 0:
                violations.append(x)

        return {
            "conjecture": "Pólya-Vermutung: L(x) ≤ 0 für alle x ≥ 2",
            "status": "WIDERLEGT (Haselgrove 1960)",
            "checked_range": f"[2, {x_max}]",
            "violations_found": len(violations),
            "first_violation": violations[0] if violations else None,
            "known_counterexample": 906316571,
            "note": "Kleinstes Gegenbeispiel ≈ 906316571 (außerhalb typischer Rechengrenze)"
        }

    def lambda_relation_to_mobius(self, n: int) -> str:
        """
        Erklärt die Beziehung zwischen λ(n) und μ(n).

        Theorem: λ(n) ist vollständig multiplikativ, μ(n) nur multiplikativ.
        Es gilt: μ(n) = Σ_{d²|n} λ(n/d²) bzw. λ(n) = Σ_{d|n} μ(d)... (Faltungsbeziehungen)

        @param n: Beispielzahl
        @return: Erklärungstext mit Werten
        @lastModified: 2026-03-12
        """
        lam_n = self.liouville(n)
        big_om = self.big_omega(n)
        return (
            f"λ({n}) = (-1)^Ω({n}) = (-1)^{big_om} = {lam_n}\n"
            f"λ ist vollständig multiplikativ: λ(mn) = λ(m)·λ(n) immer\n"
            f"μ ist nur multiplikativ: μ(mn) = μ(m)μ(n) falls gcd(m,n)=1\n"
            f"Faltung: Σ_{{d²|n}} λ(n/d²) = μ(n)"
        )
