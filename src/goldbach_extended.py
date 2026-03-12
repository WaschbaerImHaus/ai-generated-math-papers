"""
@file goldbach_extended.py
@brief Erweiterte Goldbach-Analysen: Ternäre Goldbach, Kometen, Chen, Vinogradov.
@description
    Dieses Modul implementiert erweiterte Goldbach-Analysen, insbesondere
    den ternären Satz (Helfgott 2013) und Chens Theorem.

    **Binäre Goldbach-Vermutung** (Christian Goldbach, 1742 — Conjecture, OFFEN):
        Jede gerade Zahl n ≥ 4 ist die Summe zweier Primzahlen.

    **Ternäre Goldbach-Vermutung** (Goldbach 1742 — von Helfgott 2013 BEWIESEN):
        Jede ungerade Zahl n ≥ 7 ist die Summe dreier Primzahlen.
        Beweis: Helfgott (2013) via Kreismethode + numerische Verifikation.

    **Chens Theorem** (Jingrun Chen, 1973 — BEWIESEN):
        Jede hinreichend große gerade Zahl n ist darstellbar als p + m,
        wobei p prim und m = p₁ oder m = p₁·p₂ (Produkt von ≤ 2 Primzahlen).
        Dies ist eine schwächere Version der binären Goldbach-Vermutung.

    **Goldbach-Komet** (empirisch, Conjecture):
        G(n) = #{(p,q) : p+q=n, p,q prim} für gerade n wächst "kometenartig".

    **Hardy-Littlewood Goldbach-Vermutung** (1923, Conjecture):
        G(n) ~ 2·Π₂ · ∏_{p>2, p|n} (p-1)/(p-2) · n/(log n)²
        mit Π₂ = ∏_{p>2} p(p-2)/(p-1)² ≈ 0.6601618... (Zwillingsprimkonstante)

    **Vinogradov-Theorem** (Ivan Vinogradov, 1937 — BEWIESEN):
        Jede hinreichend große ungerade Zahl ist Summe von 3 Primzahlen.
        (schwächer als Helfgott: "hinreichend groß" vs. "alle n ≥ 7")

    **Schnirelmann-Konstante** s:
        Jede positive ganze Zahl ≥ 2 ist Summe von ≤ s Primzahlen.
        Nach Oliveira e Silva (numerisch): s ≤ 4.
        Schnirelmann (1930): s ≤ 800.000 (erste endliche Schranke — BEWIESEN)

@author Michael Fuhrmann
@version 1.0
@since 2026-03-12
@lastModified 2026-03-12
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from sympy import isprime, primerange


# ===========================================================================
# HILFSFUNKTIONEN
# ===========================================================================

def _sieve(limit: int) -> np.ndarray:
    """
    Erathostenes-Sieb für schnelle Primzahlprüfung.

    @param limit: Obere Grenze (inklusiv)
    @return: Boolean-Array is_prime
    @lastModified: 2026-03-12
    """
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if is_prime[i]:
            is_prime[i * i::i] = False
    return is_prime


# ===========================================================================
# HAUPTKLASSE: GoldbachExtended
# ===========================================================================

class GoldbachExtended:
    """
    Erweiterte Goldbach-Analysen: Darstellungen, Statistiken, Heuristiken.

    Implementiert:
    - Goldbach-Darstellungen G(n) für gerade n
    - Ternäre Goldbach-Prüfung (Helfgott 2013)
    - Goldbach-Komet-Visualisierung
    - Hardy-Littlewood-Heuristik
    - Chen-Theorem-Verifikation
    - Schnirelmann-Konstante

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    # Zwillingsprimzahl-Konstante = Hardy-Littlewood Π₂ ≈ C₂
    PI2 = 0.6601618158468695739278121

    def __init__(self, limit: int = 100_000):
        """
        Initialisiert die Goldbach-Analyse mit Sieb bis limit.

        @param limit: Obere Schranke für Primzahlsieb
        @lastModified: 2026-03-12
        """
        self.limit = limit
        self._is_prime: Optional[np.ndarray] = None

    def _ensure_sieve(self) -> None:
        """Initialisiert den Primzahlsieb wenn noch nicht vorhanden."""
        if self._is_prime is None:
            self._is_prime = _sieve(self.limit)

    def is_prime(self, n: int) -> bool:
        """
        Schnelle Primzahlprüfung via Sieb.

        @param n: Zu prüfende Zahl
        @return: True wenn n prim
        @lastModified: 2026-03-12
        """
        if n > self.limit:
            return bool(isprime(n))
        self._ensure_sieve()
        return bool(self._is_prime[n])

    def goldbach_representations(self, n: int) -> List[Tuple[int, int]]:
        """
        Berechnet alle Goldbach-Darstellungen der geraden Zahl n.

        G(n) = {(p, q) : p + q = n, p ≤ q, p und q prim}

        Für ungerades n ist n = p + q nur möglich wenn p=2, dann q=n-2 muss prim sein.

        **Goldbach-Vermutung** (Conjecture, OFFEN für alle geraden n ≥ 4):
            |G(n)| ≥ 1 für alle geraden n ≥ 4.

        @param n: Zu zerlegende Zahl (idealerweise gerade, n ≥ 4)
        @return: Liste von (p, q) mit p+q=n, p ≤ q, beide prim
        @lastModified: 2026-03-12
        """
        if n < 4:
            return []
        self._ensure_sieve()
        reps = []
        # Nur p ≤ n/2 prüfen (p ≤ q)
        for p in range(2, n // 2 + 1):
            if self.is_prime(p):
                q = n - p
                if q >= p and self.is_prime(q):
                    reps.append((p, q))
        return reps

    def goldbach_count(self, n: int) -> int:
        """
        Gibt G(n) = #{(p,q) : p+q=n, p≤q prim} zurück.

        @param n: Gerade Zahl n ≥ 4
        @return: Anzahl der Goldbach-Darstellungen
        @lastModified: 2026-03-12
        """
        return len(self.goldbach_representations(n))

    def goldbach_comet(self, n_max: int = 10_000) -> List[Tuple[int, int]]:
        """
        Berechnet den Goldbach-Kometen: G(n) für alle geraden n ≤ n_max.

        Der "Komet" (Abb. der G(n)-Werte) zeigt eine charakteristische
        Struktur mit lokalen Maxima bei glatten Zahlen (vielen Primfaktoren).

        @param n_max: Maximale gerade Zahl
        @return: Liste von (n, G(n)) für alle geraden n ≥ 4
        @lastModified: 2026-03-12
        """
        comet = []
        for n in range(4, n_max + 1, 2):
            g = self.goldbach_count(n)
            comet.append((n, g))
        return comet

    def hardy_littlewood_goldbach(self, n: int) -> float:
        """
        Hardy-Littlewood-Heuristik für G(n).

        **Vermutung** (Hardy-Littlewood 1923, Conjecture A — OFFEN):
            G(n) ~ 2 · Π₂ · ∏_{p>2, p|n} (p-1)/(p-2) · n / (log n)²

        mit Π₂ = ∏_{p>2} p(p-2)/(p-1)² ≈ 0.6601618...

        @param n: Gerade Zahl n ≥ 4
        @return: Heuristische Schätzung von G(n)
        @lastModified: 2026-03-12
        """
        if n < 4 or n % 2 != 0:
            return 0.0

        singular_product = 1.0
        # Singuläres Produkt ∏_{p>2, p|n} (p-1)/(p-2)
        # Teile n durch seine ungeraden Primfaktoren
        m = n
        # Primfaktoren von n finden
        for p in range(3, min(n + 1, 10000), 2):
            if m == 1:
                break
            if m % p == 0:
                if isprime(p):
                    # Faktor (p-1)/(p-2)
                    singular_product *= (p - 1) / (p - 2)
                    while m % p == 0:
                        m //= p

        return 2.0 * self.PI2 * singular_product * n / (math.log(n) ** 2)

    def check_ternary_goldbach(self, n: int) -> Optional[Tuple[int, int, int]]:
        """
        Prüft die ternäre Goldbach-Darstellung: n = p + q + r.

        Helfgott (2013) hat die ternäre Goldbach-Vermutung vollständig bewiesen:
        Jede ungerade Zahl n ≥ 7 ist Summe von 3 Primzahlen.
        → Dies ist kein Conjecture mehr, sondern ein THEOREM.

        Für n = 5: 2+2+? nein (n-4=1 nicht prim), aber 5=2+3? Das sind nur 2.
        Für n = 7: 3+2+2 = 7 ✓

        @param n: Ungerade Zahl n ≥ 7
        @return: Tripel (p, q, r) mit p+q+r=n, alle prim, oder None
        @lastModified: 2026-03-12
        """
        if n < 5:
            return None
        self._ensure_sieve()

        # Suche p ≤ q ≤ r mit p+q+r = n
        primes_to_n = [p for p in range(2, n - 1) if self.is_prime(p)]

        for i, p in enumerate(primes_to_n):
            if p > n // 3 + 1:
                break
            remaining = n - p
            # Binäre Goldbach-Darstellung von remaining = q + r
            for q in primes_to_n[i:]:
                if q > remaining // 2 + 1:
                    break
                r = remaining - q
                if r >= q and self.is_prime(r):
                    return (p, q, r)
        return None

    def verify_ternary_goldbach_range(self, n_max: int = 1000) -> Dict:
        """
        Verifiziert die ternäre Goldbach-Vermutung für ungerade Zahlen bis n_max.

        Helfgott (2013) hat sie für ALLE ungeraden n ≥ 7 bewiesen.
        Numerische Verifikation für kleine Fälle.

        @param n_max: Maximale Grenze
        @return: Verifikationsergebnis
        @lastModified: 2026-03-12
        """
        failures = []
        count = 0
        for n in range(7, n_max + 1, 2):
            rep = self.check_ternary_goldbach(n)
            if rep is None:
                failures.append(n)
            count += 1
        return {
            "tested_up_to": n_max,
            "tested_count": count,
            "failures": failures,
            "all_verified": len(failures) == 0,
            "theorem": "Helfgott 2013: Ternäre Goldbach BEWIESEN für alle ungeraden n ≥ 7.",
        }

    def chen_theorem_check(self, n: int) -> Optional[Tuple[int, int]]:
        """
        Überprüft Chens Theorem für gerade n.

        **Chen's Theorem** (Jingrun Chen, 1973 — BEWIESEN):
            Jede hinreichend große gerade Zahl n ist darstellbar als n = p + m,
            wobei p prim und m entweder prim oder Produkt zweier Primzahlen ist.
            (Schreibweise: m ist "P₂", d.h. hat ≤ 2 Primfaktoren.)

        @param n: Gerade Zahl n ≥ 4
        @return: (p, m) mit p prim und m = prim oder Produkt zweier Primzahlen
        @lastModified: 2026-03-12
        """
        if n < 4 or n % 2 != 0:
            return None

        def is_p2(m: int) -> bool:
            """Prüft ob m ≤ 2 Primfaktoren hat (m ist P₂)."""
            if m <= 1:
                return False
            if isprime(m):
                return True
            # Prüfe ob m = p·q mit p,q prim
            if m < 4:
                return False
            for p2 in range(2, int(m ** 0.5) + 1):
                if m % p2 == 0 and isprime(p2):
                    q2 = m // p2
                    if isprime(q2):
                        return True
            return False

        for p in range(2, n):
            if self.is_prime(p):
                m = n - p
                if m >= 2 and is_p2(m):
                    return (p, m)
        return None

    def schnirelmann_density(self, n_limit: int = 100) -> Dict:
        """
        Analysiert die Schnirelmann-Dichte für Primzahlsummen.

        **Schnirelmann-Konstante s** (BEWIESEN):
            Jede positive ganze Zahl ≥ 2 ist Summe von ≤ s Primzahlen.
            Schnirelmann (1930): s < ∞ (erste endliche Schranke).
            Helfgott (2013): s ≤ 4 (für n ≥ 2, mit s=3 für ungerade n ≥ 7).

        Die Schnirelmann-Dichte σ(A) einer Menge A ⊆ ℕ ist:
            σ(A) = inf_{n≥1} #{a ∈ A : a ≤ n} / n

        @param n_limit: Bis zu welchem n geprüft wird
        @return: Dictionary mit Schnirelmann-Analyse
        @lastModified: 2026-03-12
        """
        # Prüfe minimale Primanzahl für kleine n
        prime_set = set(int(p) for p in primerange(2, n_limit + 1))
        # Füge 1 hinzu (nicht prim, aber Schnirelmann arbeitet mit Primzahlmenge)

        results = {}
        for n in range(2, min(n_limit + 1, 30)):
            if n in prime_set:
                results[n] = 1
                continue
            # Prüfe ob n = p+q (binäre Goldbach)
            found = False
            for p in prime_set:
                if p >= n:
                    break
                if (n - p) in prime_set:
                    results[n] = 2
                    found = True
                    break
            if found:
                continue
            # Prüfe ob n = p+q+r (ternäre)
            found3 = False
            for p in prime_set:
                if p >= n:
                    break
                for q in prime_set:
                    if p + q >= n:
                        break
                    r = n - p - q
                    if r in prime_set:
                        results[n] = 3
                        found3 = True
                        break
                if found3:
                    break
            if found3:
                continue
            # Dann 4 Primzahlen
            results[n] = 4  # Helfgott: ≤ 4 reicht immer

        return {
            "representation_counts": results,
            "max_primes_needed": max(results.values()) if results else 0,
            "schnirelmann_bound": 4,
            "theorem": "Schnirelmann 1930: Jede Zahl ist Summe von ≤ s < ∞ Primzahlen. "
                       "Nach Helfgott (2013): s ≤ 4.",
        }

    def goldbach_conjecture_status(self) -> Dict:
        """
        Gibt den aktuellen Status der Goldbach-Vermutungen zurück.

        @return: Status-Dictionary
        @lastModified: 2026-03-12
        """
        return {
            "binary_goldbach": {
                "statement": "Jede gerade n ≥ 4 ist Summe zweier Primzahlen",
                "status": "Conjecture — OFFEN (seit 1742)",
                "best_result": "Chens Theorem (1973): p + P₂ (fast-Goldbach)",
                "numerical_verification": "Bis 4×10^18 (Oliveira e Silva et al. 2012)",
            },
            "ternary_goldbach": {
                "statement": "Jede ungerade n ≥ 7 ist Summe dreier Primzahlen",
                "status": "BEWIESEN — Helfgott (2013)",
                "proof_method": "Kreismethode + Exponentailsummen + numerische Verifikation n < 10^30",
            },
            "chens_theorem": {
                "statement": "Jede hinreichend große gerade n = p + P₂",
                "status": "BEWIESEN — Chen (1973)",
                "method": "Große Siebbindung (large sieve) + Siebmethoden",
            },
            "vinogradov": {
                "statement": "Jede hinreichend große ungerade n ist Summe von 3 Primzahlen",
                "status": "BEWIESEN — Vinogradov (1937)",
                "bound": "Gilt für n > exp(exp(11.503)) ≈ 10^43000 (klassisch); "
                          "Helfgott reduzierte dies auf n < 8.875 × 10^30.",
            },
        }


# ===========================================================================
# KLASSE: VinogradovMethod
# ===========================================================================

class VinogradovMethod:
    """
    Implementierung der Grundideen von Vinogradovs Methode.

    Vinogradovs 3-Primzahl-Satz (1937 — BEWIESEN):
        Jede hinreichend große ungerade Zahl N ist darstellbar als
        N = p₁ + p₂ + p₃ mit p₁, p₂, p₃ prim.

    **Beweisidee (Kreismethode / Kreisverfahren)**:
    1. Darstellungsfunktion:
        r₃(N) = Σ_{p₁+p₂+p₃=N} log p₁ · log p₂ · log p₃
    2. Fourier-Darstellung via Exponentialsumme:
        r₃(N) = ∫₀¹ S(α)³ · e^{-2πiNα} dα
        mit S(α) = Σ_{p≤N} Λ(p) · e^{2πiαp}
    3. Aufteilung in Major/Minor-Bögen:
        - Major-Bögen (α nahe rationalen a/q mit kleinem q): Hauptbeitrag
        - Minor-Bögen (der Rest): Restterm → 0 via Exponentailsummen-Schranken

    @author Michael Fuhrmann
    @since 2026-03-12
    @lastModified 2026-03-12
    """

    def __init__(self):
        """
        Initialisiert die Vinogradov-Methode.
        @lastModified: 2026-03-12
        """
        pass

    def exponential_sum(self, alpha: float, N: int) -> complex:
        """
        Berechnet die Vinogradov-Exponentialsumme S(α).

        S(α) = Σ_{p≤N} log(p) · e^{2πiαp}

        Diese Summe ist das zentrale Objekt in Vinogradovs Kreismethode.
        Für α nahe einer rationalen Zahl a/q (kleine q) ist S(α) groß.
        Für α "weit von rationalen Zahlen": S(α) ist klein → Minor-Bogen-Abschätzung.

        @param alpha: Reeller Parameter α ∈ [0, 1]
        @param N: Obere Grenze der Summe
        @return: S(α) als komplexe Zahl
        @lastModified: 2026-03-12
        """
        is_p = _sieve(N)
        result = complex(0.0)
        for p in range(2, N + 1):
            if is_p[p]:
                result += math.log(p) * cmath.exp(2j * math.pi * alpha * p)
        return result

    def estimate_r3(self, N: int, num_alpha: int = 100) -> float:
        """
        Schätzt r₃(N) numerisch via Trapezregel auf [0,1].

        r₃(N) = ∫₀¹ S(α)³ · e^{-2πiNα} dα  (Parseval-Integral)

        Nur als Demonstrationsfunktion — nicht effizient für große N.

        @param N: Ungerade Zahl N (Ziel der 3-Prim-Darstellung)
        @param num_alpha: Anzahl der α-Abtastpunkte
        @return: Schätzung von r₃(N) (reeller Anteil)
        @lastModified: 2026-03-12
        """
        alphas = np.linspace(0, 1, num_alpha, endpoint=False)
        integral = complex(0.0)
        d_alpha = 1.0 / num_alpha
        for alpha in alphas:
            s = self.exponential_sum(float(alpha), N)
            integral += (s ** 3) * cmath.exp(-2j * math.pi * alpha * N) * d_alpha
        return integral.real

    def major_arc_definition(self, N: int, Q: Optional[int] = None) -> Dict:
        """
        Definiert die Major-Bögen in Vinogradovs Kreismethode.

        Major-Bögen M(a/q) = {α : |α - a/q| < N^{-2/3}} für q ≤ Q = N^{1/3}.

        Die Major-Bögen überdecken α, die nahe an rationalen a/q mit kleinem
        Nenner q liegen. Hier ist S(α) groß und liefert den Hauptbeitrag zu r₃(N).

        @param N: Obere Schranke
        @param Q: Siebbindungsparameter (Standard: N^{1/3})
        @return: Beschreibung der Major-Bögen
        @lastModified: 2026-03-12
        """
        Q_val = Q or int(N ** (1.0 / 3.0))
        arcs = []
        for q in range(1, Q_val + 1):
            for a in range(1, q + 1):
                if math.gcd(a, q) == 1:
                    center = a / q
                    width = N ** (-2.0 / 3.0)
                    arcs.append({
                        "a": a, "q": q,
                        "center": center,
                        "half_width": width,
                        "interval": (center - width, center + width),
                    })
        return {
            "Q": Q_val,
            "major_arcs": arcs[:20],  # Erste 20 zur Illustration
            "total_arcs": len(arcs),
            "explanation": (
                "Major-Bögen: α nahe a/q mit q ≤ N^{1/3} und |α-a/q| < N^{-2/3}. "
                "S(α)³ ist hier groß und liefert den Hauptterm von r₃(N)."
            ),
        }

    def minor_arc_bound_explanation(self) -> str:
        """
        Erklärt Vinogradovs Minor-Bogen-Abschätzung.

        @return: Erklärungstext
        @lastModified: 2026-03-12
        """
        return (
            "Minor-Bogen-Abschätzung (Vinogradov 1937):\n"
            "Für α auf den Minor-Bögen (außerhalb der Major-Bögen) gilt:\n"
            "    |S(α)| ≤ C · N / (ln N)^A  für jedes A > 0\n"
            "Dies folgt aus Vinogradovs Methode der Exponentialsummen\n"
            "(van der Corput / Weyl / Vinogradov-Schätzungen).\n"
            "Das Minor-Bogen-Integral ist dann:\n"
            "    |∫_{minor} S(α)³ e^{-2πiNα} dα| = o(N²/(ln N)³)\n"
            "und vernachlässigbar gegenüber dem Major-Bogen-Hauptterm:\n"
            "    ∫_{major} S(α)³ e^{-2πiNα} dα ~ S(N) · N²/(2(ln N)³)\n"
            "mit S(N) = Singuläres Serien-Produkt ≥ c > 0 für ungerade N."
        )

    def singular_series(self, N: int, prime_limit: int = 100) -> float:
        """
        Berechnet die singuläre Reihe S(N) für ungerades N.

        S(N) = ∏_{p|N} (1 − 1/(p-1)²) · ∏_{p∤N} (1 + 1/(p-1)³)

        Diese erscheint als Korrekturfaktor im Hauptterm von r₃(N):
            r₃(N) ~ S(N) · N²/(2(ln N)³)

        @param N: Ungerades N
        @param prime_limit: Obere Schranke für das Produkt
        @return: Numerische Schätzung von S(N)
        @lastModified: 2026-03-12
        """
        product = 1.0
        for p in primerange(2, prime_limit):
            p = int(p)
            if p == 2:
                continue  # Nur ungerade Primzahlen (N ist ungerade)
            if N % p == 0:
                # p teilt N: Faktor (1 − 1/(p-1)²)
                product *= 1.0 - 1.0 / ((p - 1) ** 2)
            else:
                # p teilt N nicht: Faktor (1 + 1/(p-1)³)
                product *= 1.0 + 1.0 / ((p - 1) ** 3)
        return product


# Benötigt cmath für Exponentialsumme
import cmath
