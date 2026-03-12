"""
Euler-Mascheroni-Konstante γ und Stieltjes-Konstanten – Irrationalitätstheorie.

Dieses Modul untersucht die Euler-Mascheroni-Konstante γ ≈ 0.5772156649...
und die zugehörigen Stieltjes-Konstanten γ_n. Obwohl γ mit großer Wahrscheinlichkeit
irrational ist, konnte dies bis heute nicht bewiesen werden.

Mathematischer Hintergrund:
    γ = lim_{n→∞} (Σ_{k=1}^n 1/k − ln n)
    γ_0 = γ, γ_n = lim_{N→∞} (Σ_{k=1}^N (ln k)^n / k − (ln N)^{n+1} / (n+1))
    Schranke (Papanikolaou 1997): Falls γ = p/q rational, dann q > 10^{242080}

Literatur:
    - Papanikolaou (1997): Extreme rational approximations to γ
    - Stieltjes (1885): Note sur quelques séries
    - Euler (1735): De progressionibus harmonicis observationes

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

from __future__ import annotations

import math
from typing import List, Tuple

import mpmath
import numpy as np
from sympy import EulerGamma, zeta, limit, S, log, Sum, oo, factorial, Rational


# ===========================================================================
# EULER-MASCHERONI-KONSTANTE
# ===========================================================================

class EulerMascheroniConstant:
    """
    Klasse zur hochpräzisen Berechnung und Analyse der Euler-Mascheroni-Konstante.

    Die Konstante γ tritt in der harmonischen Reihe, dem Gamma-Integral und der
    Riemann-Zeta-Funktion auf. Ihre Irrationalität ist eine der ältesten offenen
    Fragen der Zahlentheorie (Stand 2026: unbewiesen).

    Mathematische Definition:
        γ = lim_{n→∞} (H_n − ln n),  H_n = Σ_{k=1}^n 1/k

    Verbindung zu ζ(s):
        γ = lim_{s→1} (ζ(s) − 1/(s−1))

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    # Bekannter Wert (50 Stellen, als Referenz)
    GAMMA_50 = "0.57721566490153286060651209008240243104215933593992"

    def __init__(self, precision: int = 50):
        """
        Initialisiert den Rechner mit gewünschter Dezimalstellenpräzision.

        @param precision: Anzahl der Dezimalstellen (Standard: 50, max empfohlen: 100)
        """
        self.precision = precision
        mpmath.mp.dps = precision + 10  # Sicherheitspuffer

    def compute(self) -> mpmath.mpf:
        """
        Berechnet γ auf die gewünschte Präzision mittels mpmath.

        mpmath verwendet intern den Brent-McMillan-Algorithmus (O(n log n) Komplexität),
        der effizienter ist als die naive harmonische Summe.

        @return: γ als mpmath.mpf mit der eingestellten Präzision
        """
        mpmath.mp.dps = self.precision + 10
        return mpmath.euler

    def compute_via_harmonic_series(self, n_terms: int = 10000) -> float:
        """
        Berechnet γ näherungsweise über die harmonische Reihe (langsam, O(1/n)).

        Diese Methode verdeutlicht die Definition, ist aber numerisch uneffizient.
        Konvergenz: |γ − (H_n − ln n)| < 1/(2n)

        @param n_terms: Anzahl der Summanden (mehr → genauer, aber langsamer)
        @return: Näherungswert von γ als float
        """
        # Harmonische Summe H_n = 1 + 1/2 + ... + 1/n
        h_n = sum(1.0 / k for k in range(1, n_terms + 1))
        return h_n - math.log(n_terms)

    def compute_via_zeta_limit(self, epsilon: float = 1e-6) -> float:
        """
        Berechnet γ über den Grenzwert der Riemann-Zeta-Funktion.

        Mathematisch: γ = lim_{s→1} (ζ(s) − 1/(s−1))
        Numerisch: γ ≈ ζ(1+ε) − 1/ε  für kleines ε > 0

        @param epsilon: Abstand von s=1 (kleiner → genauer, aber numerisch instabiler)
        @return: Näherungswert von γ
        """
        s = 1.0 + epsilon
        zeta_val = float(mpmath.zeta(s))
        return zeta_val - 1.0 / epsilon

    def alternating_series_representation(self, n_terms: int = 200) -> float:
        """
        Alternierende Reihendarstellung für γ nach Vacca (1910).

        Formel: γ = Σ_{k=1}^∞ (−1)^k ⌊log₂ k⌋ / k
        Konvergiert langsam, aber verdeutlicht die kombinatorische Struktur.

        @param n_terms: Anzahl der Summanden
        @return: Näherungswert von γ
        """
        total = 0.0
        for k in range(1, n_terms + 1):
            # ⌊log₂ k⌋ = Bit-Länge minus 1
            floor_log2 = int(math.log2(k))
            sign = (-1) ** k
            total += sign * floor_log2 / k
        return total

    def apery_like_representation(self, n_terms: int = 100) -> float:
        """
        Apéry-artige Darstellung von γ (Euler 1735, Bailey et al.).

        Formel: γ = −Σ_{n=1}^∞ (−1)^n ln(n)/n · [Boole-Summation-Nährung]
        Benutzt: γ = Σ_{n=2}^∞ (−1)^n (⌊log₂ n⌋)/n  (Vacca-Form)

        Diese Darstellungen beweisen NICHT die Irrationalität von γ,
        zeigen aber seine tiefen Verbindungen zu harmonischen Reihen.

        @param n_terms: Anzahl der Summanden für die Näherung
        @return: Näherungswert als float
        """
        # Verwende die Formel: γ = 1 − Σ_{n=2}^∞ (ζ(n)−1)/n  (korrekte Formel)
        total = 1.0
        for n in range(2, n_terms + 1):
            # ζ(n) − 1 = Σ_{k=2}^∞ 1/k^n > 0
            zeta_n = float(mpmath.zeta(n))
            total -= (zeta_n - 1.0) / n
        return total

    def rational_approximation_bound(self) -> Tuple[str, int]:
        """
        Gibt die bekannte untere Schranke für einen rationalen Nenner zurück.

        Resultat (Papanikolaou 1997, CONJECTURE beweist Irrationalität NICHT):
            Falls γ = p/q mit ggT(p,q)=1, dann gilt: q > 10^{242080}

        Dies ist kein Beweis der Irrationalität – es schließt nur kleine Nenner aus.

        @return: Tuple (Beschreibung, Exponent) der Schranke
        """
        exponent = 242080
        description = (
            f"Papanikolaou (1997): Falls γ = p/q (rational, ggT=1), "
            f"dann gilt q > 10^{exponent}. "
            f"Dies BEWEIST NICHT die Irrationalität von γ (offene Frage, Stand 2026)."
        )
        return description, exponent

    def connection_to_gamma_function(self) -> dict:
        """
        Zeigt die Verbindung von γ zur Gamma-Funktion Γ(z).

        Formeln:
            Γ'(1) = −γ
            ln Γ(z) = −γz − ln z + Σ_{n=1}^∞ (z/n − ln(1 + z/n))
            Γ(z) = e^{−γz} / z · Π_{n=1}^∞ (1+z/n)^{−1} e^{z/n}

        @return: Dictionary mit Formel-Beschreibungen und numerischen Werten
        """
        mpmath.mp.dps = self.precision + 10
        gamma_val = mpmath.euler

        # Γ'(1) = −γ (numerisch)
        # Γ'(z) = Γ(z) · ψ(z), wobei ψ Digamma-Funktion
        digamma_1 = float(mpmath.digamma(1))  # = −γ

        return {
            "gamma": float(gamma_val),
            "digamma_at_1": digamma_1,
            "relation": "Γ'(1) = ψ(1)·Γ(1) = −γ · 1 = −γ",
            "Weierstrass_product": "Γ(z) = e^{-γz}/z · Π_{n=1}^∞ (1+z/n)^{-1} e^{z/n}",
            "verified": abs(digamma_1 + float(gamma_val)) < 1e-10
        }

    def compute_to_n_digits(self, n: int) -> str:
        """
        Berechnet γ auf genau n Dezimalstellen als String.

        @param n: Anzahl der Dezimalstellen (1 bis 100)
        @return: γ als dezimaler String der Länge n+2 (inkl. "0.")
        """
        if n < 1 or n > 100:
            raise ValueError(f"Präzision muss zwischen 1 und 100 liegen, erhalten: {n}")
        mpmath.mp.dps = n + 5
        val = mpmath.euler
        return mpmath.nstr(val, n, strip_zeros=False)

    def verify_known_digits(self) -> bool:
        """
        Überprüft die ersten 50 bekannten Stellen von γ.

        @return: True wenn Berechnung korrekt, sonst False
        """
        mpmath.mp.dps = 60
        computed = mpmath.euler
        # Vergleiche als String
        computed_str = mpmath.nstr(computed, 50, strip_zeros=False)
        # Entferne "0." Präfix für Vergleich
        known = self.GAMMA_50[2:]  # "57721566..."
        computed_digits = computed_str.replace("0.", "").replace(".", "")[:len(known)]
        return computed_digits == known


# ===========================================================================
# STIELTJES-KONSTANTEN
# ===========================================================================

class StieltjesConstants:
    """
    Berechnung und Analyse der Stieltjes-Konstanten γ_n.

    Die Stieltjes-Konstanten erscheinen in der Laurent-Entwicklung von ζ(s) um s=1:
        ζ(s) = 1/(s−1) + Σ_{n=0}^∞ (−1)^n γ_n (s−1)^n / n!

    Insbesondere gilt: γ_0 = γ (Euler-Mascheroni-Konstante)

    Mathematische Definition:
        γ_n = lim_{N→∞} (Σ_{k=1}^N (ln k)^n / k − (ln N)^{n+1} / (n+1))

    @author: Michael Fuhrmann
    @version: 1.0
    @since: 2026-03-12
    @lastModified: 2026-03-12
    """

    def __init__(self, precision: int = 30):
        """
        Initialisiert die Berechnung mit der gewünschten Präzision.

        @param precision: Dezimalstellen der Präzision (Standard: 30)
        """
        self.precision = precision
        mpmath.mp.dps = precision + 10

    def compute(self, n: int) -> mpmath.mpf:
        """
        Berechnet die n-te Stieltjes-Konstante γ_n via mpmath.

        mpmath implementiert einen effizienten Algorithmus basierend auf
        Euler-Maclaurin-Summation.

        @param n: Index der Stieltjes-Konstante (0, 1, 2, ...)
        @return: γ_n als mpmath.mpf
        """
        if n < 0:
            raise ValueError(f"Index n muss ≥ 0 sein, erhalten: {n}")
        mpmath.mp.dps = self.precision + 10
        return mpmath.stieltjes(n)

    def compute_all(self, max_n: int = 10) -> List[mpmath.mpf]:
        """
        Berechnet γ_0 bis γ_{max_n}.

        @param max_n: Maximaler Index (Standard: 10)
        @return: Liste [γ_0, γ_1, ..., γ_{max_n}]
        """
        return [self.compute(n) for n in range(max_n + 1)]

    def laurent_expansion_zeta(self, s: complex, terms: int = 6) -> complex:
        """
        Berechnet ζ(s) über die Laurent-Entwicklung (nur für s nahe 1).

        ζ(s) ≈ 1/(s−1) + Σ_{n=0}^{terms-1} (−1)^n γ_n (s−1)^n / n!

        Gültig für |s−1| < 1 (Konvergenzradius der Entwicklung).

        @param s: Komplexe Variable (s ≠ 1, |s−1| < 1 empfohlen)
        @param terms: Anzahl der Laurent-Terme (Standard: 6)
        @return: Näherungswert ζ(s)
        """
        if abs(s - 1) < 1e-10:
            raise ValueError("Laurent-Entwicklung divergiert bei s=1")

        mpmath.mp.dps = self.precision + 10
        z = mpmath.mpc(s)
        ds = z - 1  # (s − 1)

        # Polterm
        result = mpmath.mpf(1) / ds

        # Reguläre Terme
        for n in range(terms):
            gamma_n = self.compute(n)
            sign = (-1) ** n
            result += sign * gamma_n * ds**n / mpmath.factorial(n)

        return complex(result)

    def verify_gamma0_equals_euler(self) -> bool:
        """
        Überprüft, dass γ_0 = γ (Euler-Mascheroni-Konstante) gilt.

        @return: True wenn Verifikation erfolgreich
        """
        mpmath.mp.dps = 30
        gamma_0 = self.compute(0)
        euler = mpmath.euler
        return abs(float(gamma_0 - euler)) < 1e-20

    def alternating_sign_pattern(self, max_n: int = 8) -> List[int]:
        """
        Analysiert das Vorzeichenmuster der Stieltjes-Konstanten.

        Bekannt: γ_0 > 0, γ_1 < 0, höhere Terme wechseln das Vorzeichen unregelmäßig.
        Kein allgemeines Vorzeichenmuster bekannt (offene Frage).

        @param max_n: Maximaler Index
        @return: Liste der Vorzeichen (+1 oder −1)
        """
        constants = self.compute_all(max_n)
        return [1 if float(c) >= 0 else -1 for c in constants]

    def size_growth(self, max_n: int = 8) -> List[float]:
        """
        Analysiert das Wachstum |γ_n|.

        Asymptotisch: |γ_n| ~ n! / (2π)^{n+1} · ... (komplizierte Formel)
        Die Konstanten wachsen sehr schnell für große n.

        @param max_n: Maximaler Index
        @return: Liste der absoluten Werte |γ_n|
        """
        constants = self.compute_all(max_n)
        return [abs(float(c)) for c in constants]
