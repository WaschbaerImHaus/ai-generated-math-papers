"""
Agoh-Vermutung und Giuga-Vermutung: Bernoulli-Zahlen und Primzahlcharakterisierung.

Dieses Modul implementiert:
    - Bernoulli-Zahlen B_k (exakt, via sympy)
    - Agoh-Vermutung (1995): n prim ↔ n·B_{n-1} ≡ −1 (mod n)
    - Giuga-Vermutung (1950): n prim ↔ ∀ p|n: p|(n/p−1) UND (p−1)|(n/p−1)
    - Formale Äquivalenz: Agoh ↔ Giuga (via Bernoulli-Primzahl-Verteilung)
    - Giuga-Zahlen: Zusammengesetzte n mit Giuga-Eigenschaft?
    - Numerische Verifikation für n = 2..100

Mathematischer Hintergrund:
    Bernoulli-Zahlen:
        B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_4 = -1/30, ...
        Erzeugende Funktion: t/(e^t - 1) = Σ_{n≥0} B_n · t^n / n!

    Agoh-Vermutung (CONJECTURE, Agoh 1995):
        n ist prim ↔ n·B_{n-1} ≡ −1 (mod n)
        Bedeutung: Der Zähler von n·B_{n-1} + 1 ist durch n teilbar.

    Giuga-Vermutung (CONJECTURE, Giuga 1950):
        n ist prim ↔ für jeden Primteiler p von n gilt:
            (i) p | (n/p − 1)  UND  (ii) (p−1) | (n/p − 1)

    Äquivalenz (THEOREM unter Giuga-Vermutung):
        Agoh-Vermutung ↔ Giuga-Vermutung (Borwein et al. 1996)

    Giuga-Zahlen (CONJECTURE):
        Es gibt keine zusammengesetzten n mit der Giuga-Eigenschaft
        (numerisch bestätigt bis n = 10^{13800}).

    Kummer-Kongruenzen (THEOREM):
        Für gerade k, und (p-1) ∤ k: B_k/k ≡ B_{k+(p-1)}/( k+(p-1)) (mod p)

    Vorstine-Leopoldt (THEOREM):
        n ist genau dann prim, wenn n | (Σ_{k=1}^{n-1} k^{n-1} + 1)
        (→ Verbindung zu Fermat-Pseudoprimes)

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-12
@lastModified: 2026-03-12
"""

from fractions import Fraction
from typing import List, Tuple, Dict, Optional
import sympy
from sympy import bernoulli as sympy_bernoulli, isprime, factorint, Rational


def _bernoulli_fraction(n: int) -> Fraction:
    """
    Berechnet B_n als Python-Fraction (exakt).

    Nutzt sympy für zuverlässige Berechnung, konvertiert zu Fraction.

    @param n: Index der Bernoulli-Zahl n ≥ 0
    @return: B_n als exakter Bruch
    @lastModified: 2026-03-12
    """
    b = sympy_bernoulli(n)
    return Fraction(int(b.p), int(b.q))


def _gcd(a: int, b: int) -> int:
    """Euklidischer Algorithmus."""
    while b:
        a, b = b, a % b
    return abs(a)


class BernoulliNumbers:
    """
    Exakte Bernoulli-Zahlen B_0, B_1, B_2, ..., B_N via sympy.

    Bernoulli-Zahlen sind definiert durch:
        t/(e^t - 1) = Σ_{n≥0} B_n · t^n / n!

    Bekannte Werte:
        B_0 = 1
        B_1 = -1/2
        B_2 = 1/6
        B_3 = 0
        B_4 = -1/30
        B_5 = 0
        B_6 = 1/42
        B_8 = -1/30
        B_10 = 5/66

    Eigenschaften:
        - B_{2k+1} = 0 für k ≥ 1 (alle ungeraden B ≥ 3 sind 0)
        - Nenner: von-Staudt-Clausen (THEOREM): B_{2k} + Σ_{(p-1)|2k} 1/p ∈ ℤ
        - Vorzeichen: (−1)^{k+1} B_{2k} > 0

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self, max_n: int = 20):
        """
        Berechnet und cached B_0..B_max_n.

        @param max_n: Maximaler Index
        @lastModified: 2026-03-12
        """
        self.max_n = max_n
        self._cache: Dict[int, Fraction] = {}
        # Vorberechnung
        for k in range(max_n + 1):
            self._cache[k] = _bernoulli_fraction(k)

    def B(self, n: int) -> Fraction:
        """
        Gibt B_n zurück.

        @param n: Index n ≥ 0
        @return: B_n als Fraction
        @lastModified: 2026-03-12
        """
        if n < 0:
            raise ValueError(f"n={n} muss ≥ 0 sein")
        if n not in self._cache:
            self._cache[n] = _bernoulli_fraction(n)
        return self._cache[n]

    def table(self, up_to: int = 10) -> List[Tuple[int, Fraction]]:
        """
        Gibt Tabelle [(k, B_k)] für k = 0..up_to zurück.

        @param up_to: Maximaler Index
        @return: Liste von (Index, Bernoulli-Zahl)
        @lastModified: 2026-03-12
        """
        return [(k, self.B(k)) for k in range(up_to + 1)]

    def von_staudt_clausen(self, k: int) -> Fraction:
        """
        Verifiziert den Satz von von Staudt-Clausen für B_{2k}.

        THEOREM (von Staudt-Clausen, 1840):
            B_{2k} + Σ_{p: (p-1)|2k} 1/p ∈ ℤ

        Dabei summiert man über alle Primzahlen p mit (p-1) | 2k.

        @param k: Positiver Index (berechnet B_{2k})
        @return: B_{2k} + Σ 1/p (sollte ganzzahlig sein)
        @lastModified: 2026-03-12
        """
        n = 2 * k
        b_n = self.B(n)
        # Finde alle Primzahlen p mit (p-1) | 2k
        correction = Fraction(0)
        for p in range(2, n + 2):
            # Prüfe ob p prim
            if all(p % d != 0 for d in range(2, int(p**0.5) + 1)) or p == 2:
                if n % (p - 1) == 0:
                    correction += Fraction(1, p)
        result = b_n + correction
        return result


class AgohConjecture:
    """
    Agoh-Vermutung: n ist prim ↔ n·B_{n-1} ≡ −1 (mod n).

    Interpretation:
        Sei B_{n-1} = p/q (gekürzt). Dann bedeutet n·B_{n-1} ≡ -1 (mod n):
        Der Zähler von n·B_{n-1} + 1 (nach Kürzung mit q) ist durch n teilbar.

    Präzise Formulierung:
        Sei n·B_{n-1} + 1 = r/s (in niedrigsten Termen, gcd(n,s)=1 nötig).
        Die Kongruenz bedeutet: n | r.

    Agoh-Vermutung (CONJECTURE):
        Für alle n ≥ 2: n ist prim ↔ n·B_{n-1} ≡ -1 (mod n)

    Verifiziert für n ≤ 10^{13800} (kein zusammengesetztes Gegenbeispiel bekannt).

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self, bernoulli: Optional[BernoulliNumbers] = None):
        """
        Initialisiert mit optionalem BernoulliNumbers-Objekt.

        @param bernoulli: BernoulliNumbers-Instanz (wird erzeugt falls None)
        @lastModified: 2026-03-12
        """
        self._bern = bernoulli if bernoulli is not None else BernoulliNumbers(max_n=50)

    def check_agoh(self, n: int) -> bool:
        """
        Prüft die Agoh-Bedingung für n: Gilt n·B_{n-1} ≡ -1 (mod n)?

        Methode:
            1. B_{n-1} = p/q (exakter Bruch)
            2. n·B_{n-1} = n·p/q
            3. Kongruenz mod n: n·p/q + 1 ≡ 0 (mod n)?
               ⟺ n·p + q ≡ 0 (mod n·q)?
               Vereinfacht: Überprüfe ob q | gcd(q, n) = 1 vorliegt
               und n·p/q + 1 ganzzahlig ≡ 0 (mod n)

        @param n: Zu prüfende Zahl n ≥ 2
        @return: True falls n die Agoh-Bedingung erfüllt
        @lastModified: 2026-03-12
        """
        if n < 2:
            return False

        # B_{n-1} als exakter Bruch
        b = self._bern.B(n - 1)

        # n · B_{n-1} als Bruch
        nb = Fraction(n) * b  # = n * Zähler / Nenner

        # n·B_{n-1} ≡ -1 (mod n) bedeutet:
        # n·B_{n-1} + 1 ≡ 0 (mod n)
        # Sei nb = p/q (gekürzt). Dann: p/q + 1 = (p+q)/q
        # Kongruenz: (p + q)/q ≡ 0 (mod n)?
        # Nur sinnvoll falls gcd(n, q) = 1 (sonst ist Ausdruck nicht wohl-definiert mod n)

        numerator = nb.numerator + nb.denominator  # Zähler von nb + 1 = (p+q)/q
        denominator = nb.denominator

        # gcd(denominator, n) muss 1 sein für wohl-definierte Kongruenz
        g = _gcd(denominator, n)
        if g != 1:
            # Denominator hat gemeinsamen Faktor mit n
            # Kongruenz im klassischen Sinn nicht direkt anwendbar
            # Verwende: Zähler von n*B_{n-1}+1 (nach Kürzen) durch n teilbar?
            # Vollständige Kürzung:
            from math import gcd as mgcd
            g2 = mgcd(numerator, denominator)
            num_reduced = numerator // g2
            return num_reduced % n == 0

        # Standardfall: gcd(q, n) = 1 → Kongruenz p/q ≡ -1 (mod n) wohldeffiniert
        # p + q ≡ 0 (mod n) ↔ n | (p + q)
        return numerator % n == 0

    def agoh_for_range(self, n_max: int = 50) -> Dict[int, bool]:
        """
        Berechnet die Agoh-Bedingung für n = 2, ..., n_max.

        @param n_max: Obere Grenze
        @return: Dict {n: bool} (True = Agoh-Bedingung erfüllt)
        @lastModified: 2026-03-12
        """
        return {n: self.check_agoh(n) for n in range(2, n_max + 1)}

    def verify_conjecture(self, n_max: int = 50) -> Dict[str, object]:
        """
        Verifiziert die Agoh-Vermutung für alle n ≤ n_max.

        CONJECTURE: check_agoh(n) = True ↔ n ist prim

        @param n_max: Obere Grenze
        @return: Verifikationsergebnis mit Gegenbeispielen falls gefunden
        @lastModified: 2026-03-12
        """
        results = self.agoh_for_range(n_max)
        counterexamples = []

        for n, agoh_holds in results.items():
            n_is_prime = bool(isprime(n))
            if agoh_holds != n_is_prime:
                counterexamples.append({
                    "n": n,
                    "is_prime": n_is_prime,
                    "agoh_condition": agoh_holds
                })

        return {
            "conjecture": "Agoh: n prim ↔ n·B_{n-1} ≡ -1 (mod n)",
            "checked_range": f"[2, {n_max}]",
            "counterexamples": counterexamples,
            "conjecture_holds_in_range": len(counterexamples) == 0,
            "status": "CONJECTURE (kein Gegenbeispiel bis 10^{13800} bekannt)"
        }

    def agoh_bernoulli_example(self, n: int) -> Dict[str, object]:
        """
        Gibt detaillierte Berechnung der Agoh-Bedingung für eine Zahl n zurück.

        @param n: Zu prüfende Zahl
        @return: Detaillierte Berechnungsschritte
        @lastModified: 2026-03-12
        """
        b_val = self._bern.B(n - 1)
        nb = Fraction(n) * b_val
        result = self.check_agoh(n)

        return {
            "n": n,
            "n_minus_1": n - 1,
            "B_{n-1}": f"{b_val.numerator}/{b_val.denominator}",
            "n * B_{n-1}": f"{nb.numerator}/{nb.denominator}",
            "n * B_{n-1} + 1": f"{(nb + 1).numerator}/{(nb + 1).denominator}",
            "agoh_condition_holds": result,
            "is_prime": bool(isprime(n)),
            "consistent_with_conjecture": result == bool(isprime(n))
        }


class GiugaConjecture:
    """
    Giuga-Vermutung (1950): n prim ↔ ∀ Primteiler p von n gelten:
        (i) p | (n/p − 1)
        (ii) (p−1) | (n/p − 1)

    Äquivalente Formulierung (Borwein et al. 1996):
        n ist prim ↔ Σ_{k=1}^{n-1} k^{n-1} ≡ -1 (mod n)

    Giuga-Zahlen: Zusammengesetzte Zahlen mit n·B_{n-1} ≡ -1 (mod n)
    (äquivalent zur Giuga-Bedingung für alle Primteiler).

    CONJECTURE: Es gibt keine Giuga-Zahlen (verifiziert bis 10^{13800}).

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def giuga_condition_for_prime_factor(self, n: int, p: int) -> bool:
        """
        Prüft die Giuga-Bedingung für einen Primteiler p von n.

        Bedingungen:
            (i) p | (n/p − 1), d.h. (n/p − 1) ≡ 0 (mod p)
            (ii) (p−1) | (n/p − 1)

        @param n: Zusammengesetzte Zahl
        @param p: Primteiler von n
        @return: True falls beide Bedingungen erfüllt
        @lastModified: 2026-03-12
        """
        if n % p != 0:
            raise ValueError(f"{p} ist kein Teiler von {n}")
        m = n // p  # n/p
        cond1 = (m - 1) % p == 0
        cond2 = (m - 1) % (p - 1) == 0
        return cond1 and cond2

    def is_giuga_prime(self, n: int) -> bool:
        """
        Prüft, ob n die Giuga-Bedingung erfüllt (prim oder Giuga-Zahl).

        Für Primzahlen: Bedingung trivialerweise erfüllt.
        Für zusammengesetzte n: Giuga-Bedingung ⟺ n ist Giuga-Zahl (CONJECTURE: keine existieren)

        @param n: Zu prüfende Zahl n ≥ 2
        @return: True falls Giuga-Bedingung erfüllt
        @lastModified: 2026-03-12
        """
        if n < 2:
            return False
        if isprime(n):
            return True  # Alle Primzahlen erfüllen die Bedingung trivialerweise

        # Finde alle Primteiler
        prime_factors = list(factorint(n).keys())

        # Prüfe Giuga-Bedingung für jeden Primteiler
        for p in prime_factors:
            if not self.giuga_condition_for_prime_factor(n, p):
                return False
        return True

    def check_range(self, n_max: int = 100) -> Dict[str, object]:
        """
        Sucht nach Giuga-Zahlen (zusammengesetzte n mit Giuga-Eigenschaft) bis n_max.

        CONJECTURE: Keine Giuga-Zahlen existieren.

        @param n_max: Obere Grenze
        @return: Ergebnis-Dictionary
        @lastModified: 2026-03-12
        """
        giuga_composites = []
        giuga_primes = []

        for n in range(2, n_max + 1):
            if self.is_giuga_prime(n):
                if isprime(n):
                    giuga_primes.append(n)
                else:
                    giuga_composites.append(n)

        return {
            "conjecture": "Keine Giuga-Zahlen (zusammengesetzte n mit Giuga-Eigenschaft)",
            "checked_range": f"[2, {n_max}]",
            "giuga_primes_count": len(giuga_primes),
            "giuga_composites": giuga_composites,
            "conjecture_holds": len(giuga_composites) == 0,
            "status": "CONJECTURE (verifiziert bis 10^{13800})"
        }

    def giuga_sum_criterion(self, n: int) -> bool:
        """
        Prüft Giuga via Summenformel: Σ_{k=1}^{n-1} k^{n-1} ≡ -1 (mod n)?

        Äquivalente Formulierung nach Borwein et al. (1996):
            n ist prim ↔ Σ_{k=1}^{n-1} k^{n-1} ≡ -1 (mod n)

        HINWEIS: Diese Berechnung ist für große n sehr langsam (O(n²) modular).

        @param n: Zu prüfende Zahl
        @return: True falls Summen-Kriterium erfüllt
        @lastModified: 2026-03-12
        """
        # Σ_{k=1}^{n-1} k^{n-1} mod n
        total = sum(pow(k, n - 1, n) for k in range(1, n))
        return total % n == n - 1  # ≡ -1 (mod n)


class AgohGiuga:
    """
    Zusammenführung: Agoh-Vermutung ↔ Giuga-Vermutung (Äquivalenz).

    Borwein, Borwein, Borwein, Girgensohn (1996) zeigten:
        Agoh-Vermutung ↔ Giuga-Vermutung

    Beide Vermutungen sind CONJECTURES (nicht bewiesene Sätze).

    Zusammenfassung der Bernoulli-Verbindungen:
        - Wilson-Theorem (THEOREM): n prim ↔ (n-1)! ≡ -1 (mod n)
        - Fermat (THEOREM): n prim → a^{n-1} ≡ 1 (mod n) für gcd(a,n)=1
        - Giuga (CONJECTURE): n prim ↔ Σ k^{n-1} ≡ -1 (mod n)
        - Agoh (CONJECTURE): n prim ↔ n·B_{n-1} ≡ -1 (mod n)

    @author: Michael Fuhrmann
    @lastModified: 2026-03-12
    """

    def __init__(self):
        """
        Initialisiert Agoh- und Giuga-Objekte mit gemeinsamen Bernoulli-Zahlen.

        @lastModified: 2026-03-12
        """
        self.bernoulli = BernoulliNumbers(max_n=60)
        self.agoh = AgohConjecture(bernoulli=self.bernoulli)
        self.giuga = GiugaConjecture()

    def equivalence_check(self, n_max: int = 30) -> Dict[str, object]:
        """
        Überprüft die Äquivalenz Agoh ↔ Giuga für alle n ≤ n_max.

        @param n_max: Obere Grenze
        @return: Äquivalenz-Analyse
        @lastModified: 2026-03-12
        """
        divergences = []
        results = []

        for n in range(2, n_max + 1):
            agoh_val = self.agoh.check_agoh(n)
            giuga_val = self.giuga.is_giuga_prime(n)
            prime_val = bool(isprime(n))

            entry = {
                "n": n,
                "is_prime": prime_val,
                "agoh": agoh_val,
                "giuga": giuga_val,
                "consistent": agoh_val == giuga_val == prime_val
            }
            results.append(entry)

            if agoh_val != giuga_val:
                divergences.append(n)

        return {
            "conjecture": "Agoh ↔ Giuga (Borwein et al. 1996)",
            "checked_range": f"[2, {n_max}]",
            "divergences": divergences,
            "equivalence_holds_in_range": len(divergences) == 0,
            "detailed": results,
            "status": "CONJECTURE (beide equivalent aber unbewiesen)"
        }

    def bernoulli_table(self, n: int = 20) -> List[Dict[str, object]]:
        """
        Gibt eine Tabelle der Bernoulli-Zahlen B_0..B_n zurück.

        @param n: Maximaler Index
        @return: Liste von {k: B_k als String}
        @lastModified: 2026-03-12
        """
        table = []
        for k in range(n + 1):
            b = self.bernoulli.B(k)
            table.append({
                "k": k,
                "B_k": f"{b.numerator}/{b.denominator}",
                "B_k_float": float(b),
                "is_zero": b == 0,
                "numerator": b.numerator,
                "denominator": b.denominator
            })
        return table

    def comprehensive_primality_comparison(self, n: int) -> Dict[str, object]:
        """
        Vergleicht verschiedene Primzahl-Charakterisierungen für n.

        Vergleicht Wilson, Fermat (für a=2), Giuga, Agoh.

        @param n: Zu prüfende Zahl
        @return: Vergleichs-Dictionary
        @lastModified: 2026-03-12
        """
        from math import factorial

        is_p = bool(isprime(n))

        # Wilson: (n-1)! ≡ -1 (mod n)
        if n <= 20:  # Wilson für kleine n praktikabel
            wilson = factorial(n - 1) % n == n - 1
        else:
            wilson = None  # Zu aufwändig

        # Fermat für a=2: 2^{n-1} ≡ 1 (mod n) für n prim
        fermat_2 = pow(2, n - 1, n) == 1

        # Giuga
        giuga_val = self.giuga.is_giuga_prime(n)

        # Agoh
        agoh_val = self.agoh.check_agoh(n)

        # Giuga-Summe (nur für kleine n)
        giuga_sum = None
        if n <= 100:
            giuga_sum = self.giuga.giuga_sum_criterion(n)

        return {
            "n": n,
            "is_prime": is_p,
            "wilson_theorem": wilson,
            "fermat_base2": fermat_2,
            "giuga_condition": giuga_val,
            "agoh_condition": agoh_val,
            "giuga_sum_criterion": giuga_sum,
            "notes": {
                "wilson": "THEOREM (nur für p≤20 berechnet)",
                "fermat": "THEOREM (hinreichend nicht, nötig: nur Primzahlen)",
                "giuga": "CONJECTURE (Giuga 1950)",
                "agoh": "CONJECTURE (Agoh 1995)",
                "agoh_giuga_equiv": "CONJECTURE (Borwein et al. 1996)"
            }
        }
