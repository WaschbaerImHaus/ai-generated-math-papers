"""
@file normality_digits.py
@brief Analyse der Normalität von Zahlen: Ziffernverteilung, statistische Tests, Borel-Theorie.
@description
    Dieses Modul untersucht die Normalität von reellen Zahlen durch statistische
    Analyse ihrer Dezimalentwicklung.

    **Mathematische Grundlagen**:

    Eine Zahl x heißt *normal zur Basis b*, wenn jede endliche Ziffernsequenz
    der Länge k in der b-adischen Entwicklung von x mit Häufigkeit b^{-k} auftritt.

    Eine Zahl heißt *absolut normal*, wenn sie normal zu jeder Basis b≥2 ist.

    **Bekannte Ergebnisse**:
    - Champernowne-Konstante 0.123456789101112... ist BEWIESEN normal zur Basis 10
    - π, e, √2: CONJECTURE normal (extrem starke empirische Evidenz, kein Beweis!)
    - Fast alle reellen Zahlen sind normal (Borel 1909, BEWIESEN)
    - Kein bekanntes algebraisches Irrationales (außer konstruierten) ist bewiesen normal

    **Implementierte Tests**:
    - Chi²-Goodness-of-Fit: Testet Gleichverteilung der Ziffern 0–9
    - Runs-Test (Wald-Wolfowitz): Testet Unabhängigkeit aufeinanderfolgender Ziffern
    - Longest-Run-Test: Analysiert maximale Länge von Gleichlauffolgen

@author Michael Fuhrmann
@lastModified 2026-03-12
"""

import math
from typing import Optional
import mpmath
from scipy import stats
import numpy as np


# ──────────────────────────────────────────────────────────────────────────────
# Normalitätsanalyse
# ──────────────────────────────────────────────────────────────────────────────

class NormalNumberAnalysis:
    """
    @brief Statistische Analyse der Ziffernverteilung reeller Zahlen.

    @description
        Analysiert ob eine Zahl empirisch normal erscheint durch:
        1. Chi²-Test auf Gleichverteilung der Ziffern 0–9
        2. Runs-Test auf Unabhängigkeit
        3. Longest-Run-Analyse

        **WICHTIG**: Alle Tests geben nur EMPIRISCHE Evidenz.
        Der mathematische Beweis der Normalität von π, e, √2 ist OFFEN (Conjecture).

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    def __init__(self, precision: int = 1000):
        """
        @brief Initialisiert die Analyse mit gegebener Stellenanzahl.

        @param precision: Anzahl der Dezimalstellen (Standard: 1000)
        """
        self.precision = precision
        # mpmath-Präzision setzen (etwas mehr als precision für Sicherheit)
        mpmath.mp.dps = precision + 20

    def get_digits_of_pi(self, n: int) -> list:
        """
        @brief Berechnet die ersten n Dezimalziffern von π (nach dem Dezimalpunkt).

        @description
            Verwendet mpmath für hochpräzise Berechnung.
            Gibt Ziffern NACH dem Dezimalpunkt zurück (erste Stelle: 1 in 3.14...).

        @param n: Anzahl der Ziffern
        @return Liste der Ziffern als Integer 0-9
        """
        mpmath.mp.dps = n + 20
        pi_str = mpmath.nstr(mpmath.pi, n + 5, strip_zeros=False)
        # Dezimalpunkt entfernen, Vorzeichen ignorieren
        digits_str = pi_str.replace('3.', '').replace('-', '')
        return [int(d) for d in digits_str[:n] if d.isdigit()]

    def get_digits_of_e(self, n: int) -> list:
        """
        @brief Berechnet die ersten n Dezimalziffern von e (nach dem Dezimalpunkt).

        @param n: Anzahl der Ziffern
        @return Liste der Ziffern als Integer 0-9
        """
        mpmath.mp.dps = n + 20
        e_str = mpmath.nstr(mpmath.e, n + 5, strip_zeros=False)
        digits_str = e_str.replace('2.', '').replace('-', '')
        return [int(d) for d in digits_str[:n] if d.isdigit()]

    def get_digits_of_sqrt2(self, n: int) -> list:
        """
        @brief Berechnet die ersten n Dezimalziffern von √2 (nach dem Dezimalpunkt).

        @param n: Anzahl der Ziffern
        @return Liste der Ziffern als Integer 0-9
        """
        mpmath.mp.dps = n + 20
        sqrt2_str = mpmath.nstr(mpmath.sqrt(2), n + 5, strip_zeros=False)
        digits_str = sqrt2_str.replace('1.', '').replace('-', '')
        return [int(d) for d in digits_str[:n] if d.isdigit()]

    def get_digits_of_champernowne(self, n: int) -> list:
        """
        @brief Erzeugt die ersten n Dezimalziffern der Champernowne-Konstante.

        @description
            Die Champernowne-Konstante ist:
            C_10 = 0.1 2 3 4 5 6 7 8 9 10 11 12 ...

            (Aneinanderreihung aller positiven ganzen Zahlen in Dezimaldarstellung)

            **THEOREM** (Champernowne 1933, BEWIESEN):
            C_10 ist normal zur Basis 10.

            Diese Konstante dient als REFERENZ für normale Zahlen, da ihre
            Normalität bewiesen ist.

        @param n: Anzahl der benötigten Ziffern
        @return Liste der Ziffern als Integer 0-9
        """
        digits = []
        k = 1
        while len(digits) < n:
            for d in str(k):
                digits.append(int(d))
                if len(digits) >= n:
                    break
            k += 1
        return digits[:n]

    def digit_frequency(self, digits: list) -> dict:
        """
        @brief Berechnet die absolute und relative Häufigkeit jeder Ziffer 0–9.

        @param digits: Liste von Ziffern (Integer 0-9)
        @return Dictionary: {ziffer: {'count': int, 'frequency': float}}
        """
        n = len(digits)
        counts = {d: 0 for d in range(10)}

        for d in digits:
            if 0 <= d <= 9:
                counts[d] += 1

        return {
            d: {
                'count': counts[d],
                'frequency': counts[d] / n if n > 0 else 0.0
            }
            for d in range(10)
        }

    def chi_squared_test(self, digits: list) -> dict:
        """
        @brief Chi²-Goodness-of-Fit-Test auf Gleichverteilung der Ziffern 0–9.

        @description
            **Nullhypothese H_0**: Die Ziffern sind gleichverteilt auf {0,...,9}.

            Teststatistik: χ² = Σ_{d=0}^9 (O_d − E)² / E

            wobei O_d = beobachtete Häufigkeit, E = n/10 (erwartete Häufigkeit).

            Unter H_0 folgt χ² asymptotisch χ²(9) (9 Freiheitsgrade).

            **Interpretation**:
            - p-Wert > 0.05: H_0 nicht abzulehnen (keine Evidenz gegen Normalität)
            - p-Wert < 0.05: H_0 abzulehnen (Evidenz gegen Gleichverteilung)
            - p-Wert < 0.01: Starke Evidenz gegen Gleichverteilung

        @param digits: Liste von Ziffern
        @return Dictionary mit chi2-Statistik, p-Wert, Freiheitsgraden und Entscheidung
        """
        n = len(digits)
        expected = n / 10.0

        # Beobachtete Häufigkeiten
        counts = [0] * 10
        for d in digits:
            if 0 <= d <= 9:
                counts[d] += 1

        observed = np.array(counts, dtype=float)
        expected_arr = np.full(10, expected)

        # scipy chi2_contingency oder chisquare
        chi2_stat, p_value = stats.chisquare(observed, f_exp=expected_arr)

        return {
            'chi2_statistic': float(chi2_stat),
            'p_value': float(p_value),
            'degrees_of_freedom': 9,
            'n_digits': n,
            'expected_per_digit': expected,
            'observed_counts': counts,
            'reject_h0_alpha_05': p_value < 0.05,
            'reject_h0_alpha_01': p_value < 0.01,
            'interpretation': (
                'Gleichverteilung NICHT abgelehnt (p > 0.05) — konsistent mit Normalität'
                if p_value >= 0.05
                else f'Gleichverteilung abgelehnt (p = {p_value:.4f})'
            )
        }

    def runs_test(self, digits: list, threshold: Optional[int] = None) -> dict:
        """
        @brief Wald-Wolfowitz-Runs-Test auf Unabhängigkeit der Ziffernfolge.

        @description
            Der Runs-Test prüft ob aufeinanderfolgende Ziffern unabhängig sind.

            **Methode**: Binärisierung durch Vergleich mit Median (oder threshold).
            Ein "Run" ist eine maximale Folge gleicher Symbole (0 oder 1).

            **Nullhypothese H_0**: Die Ziffernfolge ist zufällig (unabhängig).

            Teststatistik: Z = (R − μ_R) / σ_R

            wobei R = Anzahl Runs, μ_R = (2n_0·n_1/(n_0+n_1)) + 1,
            σ_R² = 2n_0·n_1·(2n_0·n_1 − n_0 − n_1) / ((n_0+n_1)²·(n_0+n_1−1))

        @param digits: Liste von Ziffern
        @param threshold: Schwellenwert für Binärisierung (Standard: Median = 4)
        @return Dictionary mit Z-Statistik, p-Wert und Entscheidung
        """
        if threshold is None:
            threshold = 4  # Median der Ziffern 0–9

        # Binärisierung: 0 wenn Ziffer <= threshold, 1 sonst
        binary = [1 if d > threshold else 0 for d in digits]

        n_0 = binary.count(0)
        n_1 = binary.count(1)
        n = n_0 + n_1

        if n_0 == 0 or n_1 == 0:
            return {
                'error': 'Alle Ziffern auf einer Seite des Schwellenwerts',
                'z_statistic': None,
                'p_value': None
            }

        # Anzahl Runs zählen
        runs = 1
        for i in range(1, len(binary)):
            if binary[i] != binary[i - 1]:
                runs += 1

        # Erwartungswert und Standardabweichung der Runs
        mu_r = (2 * n_0 * n_1) / n + 1
        sigma_r_sq = (2 * n_0 * n_1 * (2 * n_0 * n_1 - n)) / (n * n * (n - 1))

        if sigma_r_sq <= 0:
            return {'error': 'Varianz nicht berechenbar', 'z_statistic': None, 'p_value': None}

        sigma_r = math.sqrt(sigma_r_sq)
        z_stat = (runs - mu_r) / sigma_r

        # Zweiseitiger p-Wert
        p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))

        return {
            'z_statistic': float(z_stat),
            'p_value': float(p_value),
            'n_runs': runs,
            'expected_runs': float(mu_r),
            'n_0': n_0,
            'n_1': n_1,
            'reject_h0_alpha_05': p_value < 0.05,
            'interpretation': (
                'Zufälligkeit NICHT abgelehnt (p > 0.05)'
                if p_value >= 0.05
                else f'Zufälligkeit abgelehnt (p = {p_value:.4f})'
            )
        }

    def longest_run_test(self, digits: list) -> dict:
        """
        @brief Analysiert die längste Folge gleicher Ziffern (Longest Run).

        @description
            Berechnet für jede Ziffer 0–9 die maximale Länge einer zusammenhängenden
            Folge dieser Ziffer.

            In einer zufälligen Folge von n Ziffern (Basis 10) erwartet man eine
            maximale Run-Länge von etwa log_{10}(n).

            **Theoretischer Erwartungswert** (approximativ):
            E[longest run of digit d] ≈ log_{10}(n)

            Längere Runs sind nicht ungewöhnlich; extrem lange würden auf
            Nicht-Normalität hindeuten.

        @param digits: Liste von Ziffern
        @return Dictionary mit maximalen Run-Längen pro Ziffer und Gesamtmaximum
        """
        n = len(digits)
        max_runs = {d: 0 for d in range(10)}

        if n == 0:
            return {'max_run_per_digit': max_runs, 'overall_max': 0, 'expected_max': 0}

        # Für jede mögliche Ziffer: längste Folge berechnen
        for target in range(10):
            current_run = 0
            max_run = 0
            for d in digits:
                if d == target:
                    current_run += 1
                    max_run = max(max_run, current_run)
                else:
                    current_run = 0
            max_runs[target] = max_run

        overall_max = max(max_runs.values())
        expected_max = math.log10(n) if n > 0 else 0

        return {
            'max_run_per_digit': max_runs,
            'overall_max': overall_max,
            'expected_max': expected_max,
            'ratio_to_expected': overall_max / expected_max if expected_max > 0 else None,
            'n_digits': n
        }

    def analyze_constant(self, name: str, n_digits: int = 1000) -> dict:
        """
        @brief Vollständige Normalitätsanalyse einer mathematischen Konstante.

        @description
            Führt alle Tests (Chi², Runs, Longest-Run) für eine Konstante durch.

            Unterstützte Konstanten:
            - 'pi': π = 3.14159...    (CONJECTURE: normal)
            - 'e': e = 2.71828...     (CONJECTURE: normal)
            - 'sqrt2': √2 = 1.41421... (CONJECTURE: normal)
            - 'champernowne': C_10     (BEWIESEN: normal zur Basis 10)

        @param name: Name der Konstante ('pi', 'e', 'sqrt2', 'champernowne')
        @param n_digits: Anzahl zu analysierender Stellen
        @return Vollständiges Analyse-Dictionary
        @raises ValueError wenn name unbekannt
        """
        # Ziffern berechnen
        if name == 'pi':
            digits = self.get_digits_of_pi(n_digits)
            status = 'CONJECTURE (nicht bewiesen)'
        elif name == 'e':
            digits = self.get_digits_of_e(n_digits)
            status = 'CONJECTURE (nicht bewiesen)'
        elif name == 'sqrt2':
            digits = self.get_digits_of_sqrt2(n_digits)
            status = 'CONJECTURE (nicht bewiesen)'
        elif name == 'champernowne':
            digits = self.get_digits_of_champernowne(n_digits)
            status = 'BEWIESEN normal zur Basis 10 (Champernowne 1933)'
        else:
            raise ValueError(
                f"Unbekannte Konstante: '{name}'. "
                f"Unterstützt: 'pi', 'e', 'sqrt2', 'champernowne'"
            )

        return {
            'constant': name,
            'normalcy_status': status,
            'n_digits_analyzed': len(digits),
            'digit_frequency': self.digit_frequency(digits),
            'chi_squared_test': self.chi_squared_test(digits),
            'runs_test': self.runs_test(digits),
            'longest_run': self.longest_run_test(digits)
        }


# ──────────────────────────────────────────────────────────────────────────────
# Normalitätsschranken (Borel-Theorie)
# ──────────────────────────────────────────────────────────────────────────────

class NormalityBounds:
    """
    @brief Theoretische Schranken und Heuristiken zur Borel-Normalitätstheorie.

    @description
        Implementiert theoretische Aspekte der Normalität:

        **Borel (1909)** [BEWIESEN]:
        Fast alle reellen Zahlen (Lebesgue-Maß 1) sind absolut normal.

        **Bekannte Normalitätsnachweise** [BEWIESEN]:
        - Champernowne-Konstante (Basis 10, 1933)
        - Copeland-Erdős-Konstante 0.23571113... (Primzahlen aneinandergereiht, 1946)
        - Stoneham-Zahlen α_{b,c} = Σ b^{-c^k·b^k} (für b,c teilerfremd, b<c)

        **Offene Probleme** [CONJECTURES]:
        - Ist π normal? (Vermutung: ja, kein Beweis)
        - Ist e normal? (Vermutung: ja, kein Beweis)
        - Ist √2 normal? (Vermutung: ja, kein Beweis)
        - Ist ln(2) normal? (Vermutung: ja, kein Beweis)
        - Gibt es algebraisch irrationale Zahlen die NICHT normal sind? (unbekannt)

    @author Michael Fuhrmann
    @lastModified 2026-03-12
    """

    def borel_measure_theorem(self) -> dict:
        """
        @brief Beschreibt Borels Maßsatz zur Normalität.

        @description
            **Satz von Borel (1909)** [BEWIESEN]:
            Die Menge der reellen Zahlen im Intervall [0,1], die NICHT normal sind,
            hat das Lebesgue-Maß Null.

            Äquivalent: Fast jede reelle Zahl (im Sinne des Maßes) ist normal.

            **Konsequenz**: "Zufällig gewählte" reelle Zahl ist fast sicher normal.

            **ABER**: Dies sagt nichts über spezifische Zahlen wie π, e, √2 aus!
            Der Beweis ist nicht konstruktiv — er zeigt Existenz normaler Zahlen,
            aber nicht, dass π etc. normal sind.

        @return Dictionary mit Theorem-Informationen
        """
        return {
            'theorem': 'Borels Maßsatz (1909)',
            'status': 'BEWIESEN',
            'statement': 'Fast alle reellen Zahlen sind absolut normal (Maß der nicht-normalen Zahlen = 0)',
            'proof_type': 'Borel-Cantelli-Lemma + Gesetz der großen Zahlen',
            'implication': 'Zufällige reelle Zahl ist fast sicher normal',
            'limitation': 'Gibt KEINE Information über spezifische Zahlen (π, e, √2)',
            'known_normal_numbers': [
                'Champernowne C_10 (Basis 10, bewiesen 1933)',
                'Copeland-Erdős-Konstante (Primzahlen, bewiesen 1946)',
                'Stoneham-Zahlen α_{b,c} (allgemein bewiesen)'
            ],
            'open_problems': [
                'π normal? (CONJECTURE)',
                'e normal? (CONJECTURE)',
                '√2 normal? (CONJECTURE)',
                'ln(2) normal? (CONJECTURE)',
                'Alle algebraischen Irrationalzahlen normal? (CONJECTURE)'
            ]
        }

    def compute_discrepancy(self, digits: list, base: int = 10) -> float:
        """
        @brief Berechnet die Diskrepanz der Ziffernverteilung von der Gleichverteilung.

        @description
            Die **Diskrepanz** misst wie stark die empirische Ziffernverteilung
            von der theoretischen Gleichverteilung abweicht:

            D_n = max_{d ∈ {0,...,b-1}} |#{i ≤ n : a_i = d}/n − 1/b|

            Eine normale Zahl hat D_n → 0 für n → ∞.

            **Berechenbare Schranken** (Conjecture für π etc.):
            Für normale Zahlen gilt D_n = O(√(log n / n)) (conjectured, nicht bewiesen für π).

        @param digits: Liste von Ziffern
        @param base: Basis (Standard: 10)
        @return Diskrepanz D_n
        """
        n = len(digits)
        if n == 0:
            return 0.0

        expected_freq = 1.0 / base
        max_deviation = 0.0

        for d in range(base):
            count = digits.count(d)
            freq = count / n
            deviation = abs(freq - expected_freq)
            max_deviation = max(max_deviation, deviation)

        return max_deviation

    def normality_bound_conjecture(self, n: int) -> dict:
        """
        @brief Berechenbare Schranken für die Diskrepanz normaler Zahlen.

        @description
            Für eine (vermutlich) normale Zahl mit n analysierten Stellen
            erwarten wir Diskrepanz D_n ~ sqrt(log(n) / n).

            **CONJECTURE** (nicht bewiesen für π, e, √2):
            Für jede normale Zahl gilt:
            D_n ≤ C · sqrt(log(n) / n) für eine Konstante C > 0.

            Für Champernowne ist diese Schranke beweisbar (explizite Konstruktion).

        @param n: Anzahl der Ziffern
        @return Dictionary mit Schranken und Heuristiken
        """
        if n <= 0:
            raise ValueError("n muss positiv sein")

        # Erwartete Diskrepanz für zufällige/normale Folge
        expected_discrepancy = math.sqrt(math.log(n) / n) if n > 1 else 1.0

        # 95%-Schranke (approximativ, aus CLT für Multinomialverteilung)
        # Standardabweichung der Häufigkeit: sqrt(p(1-p)/n) mit p=1/10
        std_deviation = math.sqrt(0.09 / n)  # sqrt(0.1 * 0.9 / n)
        upper_bound_95 = 1.96 * std_deviation  # 95% Konfidenzintervall

        return {
            'n_digits': n,
            'expected_discrepancy': expected_discrepancy,
            'theoretical_bound_95pct': upper_bound_95,
            'note': (
                'CONJECTURE für π,e,√2: D_n → 0 (normale Zahlen). '
                'Für Champernowne bewiesen.'
            ),
            'normality_criterion': 'D_n → 0 ist notwendig und hinreichend für Normalität'
        }

    def copeland_erdos_constant_digits(self, n: int) -> list:
        """
        @brief Erzeugt die ersten n Ziffern der Copeland-Erdős-Konstante.

        @description
            **Copeland-Erdős-Konstante** (BEWIESEN normal zur Basis 10, 1946):
            0.2 3 5 7 11 13 17 19 23 ...

            (Aneinanderreihung aller Primzahlen in Dezimaldarstellung)

            **Satz von Copeland und Erdős** (BEWIESEN):
            Jede Zahl der Form 0.f(1)f(2)f(3)... wobei f eine
            "hinreichend dichte" Folge natürlicher Zahlen ist, ist normal zur Basis 10.
            Primzahlen erfüllen diese Dichtebedingung (Primzahlsatz).

        @param n: Anzahl der benötigten Ziffern
        @return Liste der Ziffern als Integer 0-9
        """
        digits = []
        p = 2  # Startprime

        while len(digits) < n:
            for d in str(p):
                digits.append(int(d))
                if len(digits) >= n:
                    break
            # Nächste Primzahl
            p = self._next_prime(p)

        return digits[:n]

    def _next_prime(self, p: int) -> int:
        """
        @brief Berechnet die nächste Primzahl nach p.

        @param p: Aktuelle Primzahl
        @return Nächste Primzahl größer als p
        """
        candidate = p + 1
        while True:
            if all(candidate % i != 0 for i in range(2, int(math.isqrt(candidate)) + 1)):
                return candidate
            candidate += 1
