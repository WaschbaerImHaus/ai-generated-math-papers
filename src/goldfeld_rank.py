"""
Goldfeld-Rangvermutung und Selmer-Gruppen-Analyse elliptischer Kurven.

Dieses Modul implementiert:
    - Goldfeld-Rangvermutung: Durchschnittlicher Rang = 1/2 (vermutet 1979)
    - Bhargava-Shankar-Schranken (2013/2015): avg. Rang ≤ 7/6 (2-Selmer)
    - Rang-Verteilung für Kurven y² = x³ + ax + b
    - Naiver 2-Selmer-Rang über quadratische Punkte
    - Verbindung zu Birch & Swinnerton-Dyer

Mathematischer Hintergrund:
    Elliptische Kurven E/ℚ: y² = x³ + Ax + B, Diskriminante Δ ≠ 0

    Birch & Swinnerton-Dyer (CONJECTURE):
        rang(E(ℚ)) = ord_{s=1} L(E, s)

    Goldfeld-Vermutung (CONJECTURE, 1979):
        (1/N) Σ_{|a|,|b|≤N} rang(E_{a,b}) → 1/2 für N → ∞

    Bhargava-Shankar (THEOREM, 2013):
        avg. 2-Selmer-Rang ≤ 3  →  avg. Rang ≤ 3/2
    Bhargava-Shankar (THEOREM, 2015):
        avg. 3-Selmer-Rang ≤ 2  →  avg. Rang ≤ 7/6 ≈ 1.167

    Heegner-Punkte (THEOREM, Gross-Zagier + Kolyvagin):
        Falls L(E,1) ≠ 0 → rang = 0
        Falls L'(E,1) ≠ 0, ord=1 → rang = 1

    Goldfeld (1983): h(-D) ≥ (1/7000) log D für Heegner-Punkte verwendet.
    h(-163) = 1: einziger imaginär-quadratischer Körper der Klassenzahl 1
    mit |Diskriminante| = 163 (Baker-Heegner-Stark, 1966–1967).

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

import math
from typing import List, Tuple, Dict, Optional
from fractions import Fraction


def _gcd(a: int, b: int) -> int:
    """Euklidischer Algorithmus für ggT."""
    while b:
        a, b = b, a % b
    return abs(a)


def _discriminant(a: int, b: int) -> int:
    """
    Berechnet die Diskriminante Δ = -16(4A³ + 27B²) der Weierstraß-Kurve y² = x³ + Ax + B.

    Δ ≠ 0: Kurve ist nicht-singulär (echte elliptische Kurve).

    @param a: Koeffizient A
    @param b: Koeffizient B
    @return: Diskriminante Δ
    @lastModified: 2026-03-12
    """
    return -16 * (4 * a ** 3 + 27 * b ** 2)


def _is_quadratic_residue(a: int, p: int) -> bool:
    """
    Prüft, ob a ein quadratischer Rest mod p ist (Euler-Kriterium).

    a ist QR mod p ↔ a^{(p-1)/2} ≡ 1 (mod p)

    @param a: Ganzzahl
    @param p: Ungerade Primzahl
    @return: True wenn a ein quadratischer Rest mod p ist
    @lastModified: 2026-03-12
    """
    if p == 2:
        return True
    r = pow(a % p, (p - 1) // 2, p)
    return r == 1


class EllipticCurveRank:
    """
    Repräsentiert eine elliptische Kurve E: y² = x³ + Ax + B über ℚ.

    Berechnet grundlegende Eigenschaften und eine heuristische Rang-Abschätzung.

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self, a: int, b: int):
        """
        Initialisiert die elliptische Kurve y² = x³ + ax + b.

        @param a: Koeffizient A in der Weierstraß-Form
        @param b: Koeffizient B in der Weierstraß-Form
        @raises ValueError: Falls die Diskriminante 0 ist (singuläre Kurve)
        @lastModified: 2026-03-12
        """
        self.a = a
        self.b = b
        self.disc = _discriminant(a, b)
        if self.disc == 0:
            raise ValueError(f"Singuläre Kurve: Δ=0 für a={a}, b={b}")

    def __repr__(self) -> str:
        return f"EllipticCurve(y²=x³+{self.a}x+{self.b}, Δ={self.disc})"

    def is_on_curve(self, x: int, y: int) -> bool:
        """
        Prüft ob der Punkt (x, y) ∈ E(ℚ) liegt.

        @param x: x-Koordinate
        @param y: y-Koordinate
        @return: True wenn y² = x³ + ax + b
        @lastModified: 2026-03-12
        """
        return y * y == x * x * x + self.a * x + self.b

    def points_mod_p(self, p: int) -> List[Tuple[int, int]]:
        """
        Findet alle affinen Punkte auf E(𝔽_p) = E mod p.

        Brute-Force über alle x, y ∈ 𝔽_p.

        @param p: Primzahl p (der endliche Körper 𝔽_p)
        @return: Liste aller affinen Punkte (x, y) mod p
        @lastModified: 2026-03-12
        """
        pts = []
        for x in range(p):
            rhs = (x ** 3 + self.a * x + self.b) % p
            for y in range(p):
                if (y * y) % p == rhs:
                    pts.append((x, y))
        return pts

    def ap(self, p: int) -> int:
        """
        Berechnet a_p = p + 1 − #E(𝔽_p) (Frobenius-Spur).

        Hasse-Schranke (THEOREM): |a_p| ≤ 2√p

        @param p: Primzahl
        @return: Frobenius-Spur a_p
        @lastModified: 2026-03-12
        """
        # Zähle affine Punkte + unendlichen Punkt
        count = len(self.points_mod_p(p)) + 1
        return p + 1 - count

    def naive_rank_bound(self, prime_limit: int = 50) -> int:
        """
        Heuristik: Schätzt Rang via 2-Descent-Näherung ab.

        HINWEIS: Dies ist nur eine grobe Heuristik und keine exakte Berechnung.
        Exakte Rang-Bestimmung erfordert vollständigen 2-Descent oder L-Funktionen.

        Methode: Zähle Primzahlen p ≤ prime_limit, für die mehr Punkte als
        erwartet vorhanden sind (Indiz für positiven Rang).

        @param prime_limit: Obergrenze für Primzahlen
        @return: Heuristischer Rang-Schätzer ∈ {0, 1, 2, ...}
        @lastModified: 2026-03-12
        """
        primes = [p for p in range(2, prime_limit + 1) if all(p % d != 0 for d in range(2, int(p**0.5) + 1))]
        # Einfache Heuristik: Wenn a_p häufig klein ist → höherer Rang
        # BSW-Heuristik: Π_{p≤N} #E(𝔽_p)/p ~ c·(log N)^r, r = Rang
        product = 1.0
        for p in primes:
            if self.disc % p == 0:
                continue  # schlechte Reduktion
            n_p = len(self.points_mod_p(p)) + 1  # + Punkt im Unendlichen
            product *= n_p / p

        if product < 1:
            return 0
        # Grobe Schätzung: r ≈ 2·log(product)/log(log(prime_limit))
        log_prod = math.log(product) if product > 0 else 0
        r = max(0, round(2 * log_prod / math.log(math.log(prime_limit + 2) + 1)))
        return r

    def two_selmer_rank_bound(self) -> int:
        """
        Einfache obere Schranke für den 2-Selmer-Rang.

        Die 2-Selmer-Gruppe S^{(2)}(E/ℚ) enthält E(ℚ)/2E(ℚ) als Untergruppe.
        2-Selmer-Rang ≥ Rang(E) (obere Schranke für Rang).

        Methode: Zähle Quadratreste der rechten Seite bei kleinen x-Werten.
        Gibt nur eine grobe Schranke, keinen exakten Wert.

        @return: Obere Schranke für Rang via 2-Selmer
        @lastModified: 2026-03-12
        """
        # Zähle rationale Punkte für kleine |x|
        count = 0
        for x_num in range(-20, 21):
            for x_den in range(1, 5):
                x = Fraction(x_num, x_den)
                rhs = x**3 + self.a * x + self.b
                if rhs > 0:
                    sqrt_rhs = math.sqrt(float(rhs))
                    # Prüfe ob rhs ein vollständiges Quadrat einer rationalen Zahl ist
                    # Einfache Ganzzahlprüfung für kleine Werte
                    if x_den == 1:
                        y_int = round(sqrt_rhs)
                        if abs(y_int**2 - float(rhs)) < 1e-6 and y_int > 0:
                            count += 1
                elif rhs == 0:
                    count += 1  # x-Achsen-Punkt

        # Sehr grobe Schranke: count/2 als Rang-Indikator
        return min(count // 2, 10)


class GoldfeldRankConjecture:
    """
    Analysiert die Goldfeld-Rangvermutung: Durchschnittlicher Rang = 1/2.

    CONJECTURE (Goldfeld 1979):
        Für die Familie {E_{a,b}: y² = x³ + ax + b, a,b ∈ ℤ} gilt:
        lim_{N→∞} (1/N²) Σ_{|a|,|b|≤N} rang(E_{a,b}) = 1/2

    Bekannte Schranken (THEOREME):
        Bhargava-Shankar (2013): avg. 2-Selmer-Rang ≤ 3, daher avg. Rang ≤ 3/2
        Bhargava-Shankar (2015): avg. 3-Selmer-Rang ≤ 2, daher avg. Rang ≤ 7/6
        Bhargava-Shankar (2015): avg. 5-Selmer → avg. Rang ≤ 1.17

    Heegner-Punkte / Gross-Zagier (1986) (THEOREM):
        Wenn ord_{s=1} L(E,s) = 1, dann rang(E(ℚ)) = 1 und Sha(E/ℚ) endlich.

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    # Bekannte Bhargava-Shankar-Schranken als Klassenattribute
    BHARGAVA_SHANKAR_2SELMER = Fraction(3, 2)    # 2013: avg. Rang ≤ 3/2
    BHARGAVA_SHANKAR_3SELMER = Fraction(7, 6)    # 2015: avg. Rang ≤ 7/6
    GOLDFELD_CONJECTURE = Fraction(1, 2)          # Vermuteter Grenzwert

    def __init__(self, a_range: int = 10, b_range: int = 10):
        """
        Initialisiert die Analyse mit einem Bereich von Kurven.

        @param a_range: Koeffizienten a ∈ [-a_range, a_range]
        @param b_range: Koeffizienten b ∈ [-b_range, b_range]
        @lastModified: 2026-03-12
        """
        self.a_range = a_range
        self.b_range = b_range
        self._curves: List[EllipticCurveRank] = []
        self._rank_distribution: Dict[int, int] = {}

    def build_curve_family(self) -> List[EllipticCurveRank]:
        """
        Erstellt die Familie der elliptischen Kurven für a,b ∈ [-range, range].

        Singuläre Kurven (Δ=0) werden übersprungen.

        @return: Liste nicht-singulärer elliptischer Kurven
        @lastModified: 2026-03-12
        """
        self._curves = []
        for a in range(-self.a_range, self.a_range + 1):
            for b in range(-self.b_range, self.b_range + 1):
                if _discriminant(a, b) != 0:
                    self._curves.append(EllipticCurveRank(a, b))
        return self._curves

    def rank_distribution(self) -> Dict[int, int]:
        """
        Berechnet die Rang-Verteilung für die Kurven-Familie.

        Nutzt naive_rank_bound() als Heuristik.

        @return: Dict {Rang: Anzahl Kurven}
        @lastModified: 2026-03-12
        """
        if not self._curves:
            self.build_curve_family()

        self._rank_distribution = {}
        for curve in self._curves:
            r = curve.naive_rank_bound()
            self._rank_distribution[r] = self._rank_distribution.get(r, 0) + 1
        return self._rank_distribution

    def average_rank_estimate(self) -> float:
        """
        Schätzt den Durchschnittsrang der Familie.

        CONJECTURE (Goldfeld): Sollte → 1/2 konvergieren.
        Bhargava-Shankar: Liegt ≤ 7/6 ≈ 1.167.

        @return: Heuristischer Durchschnittsrang
        @lastModified: 2026-03-12
        """
        if not self._rank_distribution:
            self.rank_distribution()

        total = sum(self._rank_distribution.values())
        weighted = sum(r * cnt for r, cnt in self._rank_distribution.items())
        if total == 0:
            return 0.0
        return weighted / total

    def bhargava_shankar_comparison(self) -> Dict[str, object]:
        """
        Vergleicht die berechneten Schranken mit Bhargava-Shankar-Resultaten.

        @return: Vergleichs-Dictionary
        @lastModified: 2026-03-12
        """
        avg = self.average_rank_estimate()
        return {
            "family": f"y²=x³+ax+b, a,b∈[-{self.a_range},{self.a_range}]×[-{self.b_range},{self.b_range}]",
            "n_curves": len(self._curves) if self._curves else "nicht berechnet",
            "average_rank_heuristic": avg,
            "goldfeld_conjecture": float(self.GOLDFELD_CONJECTURE),
            "bhargava_shankar_2selmer_bound": float(self.BHARGAVA_SHANKAR_2SELMER),
            "bhargava_shankar_3selmer_bound": float(self.BHARGAVA_SHANKAR_3SELMER),
            "rank_distribution": self._rank_distribution.copy(),
            "notes": [
                "CONJECTURE (Goldfeld 1979): avg. Rang → 1/2",
                "THEOREM (Bhargava-Shankar 2013): avg. Rang ≤ 3/2",
                "THEOREM (Bhargava-Shankar 2015): avg. Rang ≤ 7/6",
                "Heuristik hier unzuverlässig für kleine |a|,|b|"
            ]
        }

    def heegner_point_example(self) -> Dict[str, object]:
        """
        Erklärt Heegner-Punkte und den Zusammenhang mit der Klassenzahl h(-163) = 1.

        Gross-Zagier-Formel (THEOREM, 1986):
            Falls L(E, 1) = 0 und L'(E, 1) ≠ 0:
            h(y_K)² = (||f||² / u²) · L'(E, 1) · √|D| / (2π)^2

        Goldfeld (1983) nutzte Modulkurven und Heegner-Punkte, um eine effektiv
        berechenbare untere Schranke für h(-D) zu beweisen:
            h(-D) ≥ (1/7000) · log D

        h(-163) = 1: Einziger imaginär-quadratischer Körper Q(√-163) der Klassenzahl 1
        außerhalb der 9 Heegner-Zahlen. (Baker-Heegner-Stark, 1966–67)

        @return: Erklärungstext als Dictionary
        @lastModified: 2026-03-12
        """
        return {
            "heegner_number_163": 163,
            "class_number_h_minus_163": 1,
            "explanation": (
                "h(-163) = 1: Der imaginär-quadratische Körper Q(√-163) "
                "hat Klassenzahl 1 (eindeutige Primfaktorzerlegung). "
                "Dies ist die größte der 9 Heegner-Zahlen "
                "{1,2,3,7,11,19,43,67,163} (Baker-Heegner-Stark 1966-67)."
            ),
            "goldfeld_lower_bound": "h(-D) ≥ (1/7000)·log(D) [THEOREM, Goldfeld 1983]",
            "gross_zagier": "Rang 1 ↔ L'(E,1) ≠ 0 [THEOREM, Gross-Zagier 1986]",
            "kolyvagin": "Rang 0,1 ↔ Sha(E/Q) endlich [THEOREM, Kolyvagin 1988]"
        }


class SelmerGroupAnalysis:
    """
    Analyse der 2-Selmer-Gruppe einer elliptischen Kurve.

    Die 2-Selmer-Gruppe S^{(2)}(E/ℚ) passt in die exakte Sequenz:
        0 → E(ℚ)/2E(ℚ) → S^{(2)}(E/ℚ) → Ш(E/ℚ)[2] → 0

    wobei Ш(E/ℚ) die Tate-Shafarevich-Gruppe ist.

    Rang(E(ℚ)) ≤ dim S^{(2)}(E/ℚ) - (Beitrag von Ш)

    CONJECTURE (BSD): Ш(E/ℚ) ist endlich für alle E.

    Bhargava-Shankar (2013) zeigten:
        avg. dim S^{(2)} = 3  (gemittelt über alle elliptischen Kurven geordnet nach Höhe)

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self, curve: EllipticCurveRank):
        """
        Initialisiert die Selmer-Analyse für eine gegebene Kurve.

        @param curve: EllipticCurveRank-Objekt
        @lastModified: 2026-03-12
        """
        self.curve = curve

    def two_descent_primes(self, bound: int = 100) -> List[int]:
        """
        Identifiziert relevante Primzahlen für den 2-Descent.

        Für den 2-Descent sind relevant:
        - p = 2
        - Primzahlen p | Δ (schlechte Reduktion)
        - Primzahlen p | (4A³ + 27B²)

        @param bound: Obere Schranke für Primzahlen
        @return: Liste relevanter Primzahlen
        @lastModified: 2026-03-12
        """
        relevant = {2}
        # Primteiler der Diskriminante
        d = abs(self.curve.disc)
        k = d
        p = 2
        while p * p <= k:
            if k % p == 0:
                relevant.add(p)
                while k % p == 0:
                    k //= p
            p += 1
        if k > 1:
            relevant.add(k)

        return sorted(p for p in relevant if p <= bound)

    def selmer_rank_upper_bound(self) -> int:
        """
        Berechnet eine obere Schranke für dim S^{(2)}(E/ℚ).

        Methode: 2-Descent via Primteiler-Analyse.
        Die genaue Berechnung erfordert lokale Lösbarkeit und Galois-Kohomologie.
        Diese Implementierung gibt eine einfache obere Schranke.

        THEOREM (Bhargava-Shankar):
            E[avg. dim S^{(2)}] = 3 (geordnet nach Höhe H(E) = max(4|A|³, 27B²))

        @return: Obere Schranke für Selmer-Rang
        @lastModified: 2026-03-12
        """
        # Anzahl der 2-Descent-Primes als Indikator
        primes = self.two_descent_primes()
        # Selmer-Rang ≤ |schlechte Primes| + 2 (grobe Schranke)
        return len(primes) + 2

    def torsion_subgroup_points(self, bound: int = 20) -> List[Tuple[int, int]]:
        """
        Sucht Torsionspunkte der Kurve mit kleinen ganzzahligen Koordinaten.

        Nagell-Lutz-Theorem (THEOREM):
            Jeder ganzzahlige Torsionspunkt (x,y) mit y ≠ 0 erfüllt y² | Δ.

        @param bound: Suchbereich für x,y-Koordinaten
        @return: Liste möglicher Torsionspunkte
        @lastModified: 2026-03-12
        """
        torsion = []
        disc_abs = abs(self.curve.disc)
        # Kandidaten: y² | Δ
        y_candidates = [0] + [y for y in range(1, bound + 1)
                               if disc_abs % (y * y) == 0 or y * y <= bound * 2]

        for x in range(-bound, bound + 1):
            rhs = x**3 + self.curve.a * x + self.curve.b
            if rhs < 0:
                continue
            for y in y_candidates:
                if y * y == rhs:
                    torsion.append((x, y))
                    if y > 0:
                        torsion.append((x, -y))
        return torsion

    def summary(self) -> Dict[str, object]:
        """
        Gibt eine Zusammenfassung der Selmer-Analyse zurück.

        @return: Zusammenfassungs-Dictionary
        @lastModified: 2026-03-12
        """
        return {
            "curve": str(self.curve),
            "discriminant": self.curve.disc,
            "descent_primes": self.two_descent_primes(),
            "selmer_rank_bound": self.selmer_rank_upper_bound(),
            "torsion_points": self.torsion_subgroup_points(),
            "bsd_connection": (
                "BSD CONJECTURE: rang(E(ℚ)) = ord_{s=1} L(E,s). "
                "Dim S^(2) ≥ rang ≥ 0."
            ),
            "bhargava_shankar": (
                "THEOREM: avg. dim S^(2) = 3 über alle E/Q "
                "(Bhargava-Shankar 2013)"
            )
        }
