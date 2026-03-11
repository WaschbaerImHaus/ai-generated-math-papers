"""
@file statistics_math.py
@brief Statistik-Modul für das specialist-maths Projekt.
@description
    Enthält deskriptive Statistik, Wahrscheinlichkeitsverteilungen,
    Hypothesentests, Bayes-Theorem und Monte-Carlo-Simulation.

    Implementierte Funktionen:
    - Deskriptive Statistik: mean, median, mode, variance, std_dev,
      quartiles, iqr, skewness, kurtosis
    - Normalverteilung: PDF, CDF, PPF (Quantilfunktion)
    - Binomialverteilung: PMF, CDF
    - Poisson-Verteilung: PMF, CDF
    - Hypothesentests: t-Test (ein/zwei Stichproben), Chi-Quadrat-Test
    - Bayes-Theorem
    - Monte-Carlo-Simulation: Pi-Schätzung

@author Michael Fuhrmann
@date 2026-03-05
@lastModified 2026-03-10
"""

import math
import random
import numpy as np
from typing import Any


# ─────────────────────────────────────────────────────────────────────────────
# Hilfsfunktionen
# ─────────────────────────────────────────────────────────────────────────────

def _erf(x: float) -> float:
    """
    Berechnet die Gauss'sche Fehlerfunktion erf(x) via Horner-Approximation.

    Die Fehlerfunktion ist definiert als:
      erf(x) = (2/sqrt(pi)) * integral_0^x exp(-t^2) dt

    Approximation nach Abramowitz & Stegun, Formel 7.1.26.

    @param x Eingabewert
    @return Näherungswert von erf(x), Fehler < 1.5e-7
    @lastModified 2026-03-10
    """
    # Vorzeichenbehandlung: erf ist ungerade
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # Koeffizienten der Approximation (Abramowitz & Stegun)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (((((1.061405429 * t - 1.453152027) * t)
                  + 1.421413741) * t - 0.284496736) * t
                + 0.254829592) * t * math.exp(-x * x)
    return sign * y


def _erfinv(p: float) -> float:
    """
    Berechnet die inverse Fehlerfunktion erfinv(p).

    Verwendet eine Rational-Approximation (Peter Acklam's Methode).

    @param p Eingabewert im Bereich (-1, 1)
    @return x mit erf(x) = p
    @lastModified 2026-03-10
    """
    # Umwandlung: erfinv(p) = ppf((1+p)/2) / sqrt(2) für Normalverteilung
    # Wir verwenden die Probit-Funktion über rational approximation
    return _probit((1.0 + p) / 2.0) / math.sqrt(2.0)


def _probit(p: float) -> float:
    """
    Berechnet die Probit-Funktion (Quantilfunktion der Standardnormalverteilung).

    Rational-Approximation nach Peter Acklam (max. Fehler: 1.15e-9).

    @param p Wahrscheinlichkeit im Bereich (0, 1)
    @return z-Wert mit CDF(z) = p
    @lastModified 2026-03-10
    """
    if p <= 0 or p >= 1:
        raise ValueError(f"p muss im Bereich (0,1) liegen, nicht {p}")

    # Koeffizienten der rationalen Approximation
    a = [-3.969683028665376e+01,  2.209460984245205e+02,
         -2.759285104469687e+02,  1.383577518672690e+02,
         -3.066479806614716e+01,  2.506628277459239e+00]
    b = [-5.447609879822406e+01,  1.615858368580409e+02,
         -1.556989798598866e+02,  6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
          4.374664141464968e+00,  2.938163982698783e+00]
    d = [7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00]

    p_low  = 0.02425
    p_high = 1 - p_low

    if p_low <= p <= p_high:
        # Mittlerer Bereich: rationale Approximation
        q = p - 0.5
        r = q * q
        num = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q
        den = (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)
        return num / den
    elif p < p_low:
        # Linker Rand
        q = math.sqrt(-2 * math.log(p))
        num = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5])
        den = ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
        return num / den
    else:
        # Rechter Rand: Symmetrie
        return -_probit(1 - p)


def _gamma_inc_reg(a: float, x: float) -> float:
    """
    Berechnet die regularisierte untere unvollständige Gammafunktion P(a,x).

    Verwendet Reihenentwicklung für x < a+1, Kettenbruch für x >= a+1.

    @param a Formparameter (a > 0)
    @param x Obere Grenze (x >= 0)
    @return P(a,x) = gamma(a,x) / Gamma(a)
    @lastModified 2026-03-10
    """
    if x < 0:
        raise ValueError("x muss >= 0 sein")
    if x == 0:
        return 0.0

    # Log-Gammafunktion via Lanczos-Approximation
    log_gamma_a = math.lgamma(a)

    if x < a + 1:
        # Reihenentwicklung
        ap = a
        delta = 1.0 / a
        total = delta
        for _ in range(200):
            ap += 1
            delta *= x / ap
            total += delta
            if abs(delta) < abs(total) * 1e-12:
                break
        return total * math.exp(-x + a * math.log(x) - log_gamma_a)
    else:
        # Kettenbruch (Lentz-Algorithmus)
        b = x + 1 - a
        c = 1.0 / 1e-30
        d = 1.0 / b
        h = d
        for i in range(1, 201):
            an = -i * (i - a)
            b += 2
            d = an * d + b
            if abs(d) < 1e-30:
                d = 1e-30
            c = b + an / c
            if abs(c) < 1e-30:
                c = 1e-30
            d = 1.0 / d
            delta = d * c
            h *= delta
            if abs(delta - 1) < 1e-12:
                break
        # Obere unvollständige -> 1 - P(a,x)
        return 1.0 - math.exp(-x + a * math.log(x) - log_gamma_a) * h


def _t_cdf(t: float, df: int) -> float:
    """
    Berechnet die CDF der t-Verteilung P(T <= t) mit df Freiheitsgraden.

    Verwendet die Beziehung zur regularisierten Betafunktion.

    @param t t-Statistik
    @param df Freiheitsgrade
    @return P(T <= t)
    @lastModified 2026-03-10
    """
    x = df / (df + t * t)
    # Regularisierte Betafunktion I_x(df/2, 0.5)
    p = _beta_inc_reg(df / 2.0, 0.5, x)
    if t >= 0:
        return 1.0 - p / 2.0
    else:
        return p / 2.0


def _beta_inc_reg(a: float, b: float, x: float) -> float:
    """
    Berechnet die regularisierte unvollständige Betafunktion I_x(a,b).

    Verwendet den Algorithmus von Numerical Recipes (Kettenbruch).

    @param a Erster Formparameter
    @param b Zweiter Formparameter
    @param x Obere Grenze (0 <= x <= 1)
    @return I_x(a,b)
    @lastModified 2026-03-10
    """
    if x < 0 or x > 1:
        raise ValueError(f"x muss in [0,1] liegen, nicht {x}")
    if x == 0:
        return 0.0
    if x == 1:
        return 1.0

    # Log des Beta-Koeffizienten
    log_beta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1 - x) - log_beta)

    # Für bessere Konvergenz: Symmetrierelation I_x(a,b) = 1 - I_{1-x}(b,a)
    if x > (a + 1) / (a + b + 2):
        return 1.0 - _beta_inc_reg(b, a, 1 - x)

    # Lentz-Kettenbruch
    def cont_frac() -> float:
        qab = a + b
        qap = a + 1
        qam = a - 1
        c = 1.0
        d = 1.0 - qab * x / qap
        if abs(d) < 1e-30:
            d = 1e-30
        d = 1.0 / d
        h = d
        for m in range(1, 101):
            m2 = 2 * m
            # Gerades m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1 + aa * d
            if abs(d) < 1e-30:
                d = 1e-30
            c = 1 + aa / c
            if abs(c) < 1e-30:
                c = 1e-30
            d = 1.0 / d
            h *= d * c
            # Ungerades m
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1 + aa * d
            if abs(d) < 1e-30:
                d = 1e-30
            c = 1 + aa / c
            if abs(c) < 1e-30:
                c = 1e-30
            d = 1.0 / d
            delta = d * c
            h *= delta
            if abs(delta - 1) < 1e-12:
                break
        return h

    return front * cont_frac() / a


def _chi2_cdf(x: float, df: int) -> float:
    """
    Berechnet die CDF der Chi-Quadrat-Verteilung P(X <= x) mit df Freiheitsgraden.

    Verwendet die regularisierte Gammafunktion: CDF = P(df/2, x/2).

    @param x Eingabewert (x >= 0)
    @param df Freiheitsgrade
    @return P(X <= x)
    @lastModified 2026-03-10
    """
    if x < 0:
        return 0.0
    return _gamma_inc_reg(df / 2.0, x / 2.0)


def _normal_cdf_approx(x: float) -> float:
    """
    Berechnet eine schnelle Approximation der Standard-Normalverteilungs-CDF.

    Hilfsfunktion für interne Berechnungen, wo Geschwindigkeit wichtiger
    als höchste Genauigkeit ist.

    @param x Standardisierter z-Wert
    @return Näherungswert von P(Z <= x) für Z ~ N(0,1)
    @lastModified 2026-03-10
    """
    return 0.5 * (1.0 + _erf(x / math.sqrt(2.0)))


# ─────────────────────────────────────────────────────────────────────────────
# Deskriptive Statistik
# ─────────────────────────────────────────────────────────────────────────────

def mean(data: list[float]) -> float:
    """
    Berechnet den arithmetischen Mittelwert einer Datenmenge.

    Formel: mu = (1/n) * sum(x_i)

    Implementierung mit np.mean() für optimale Performance bei großen Datensätzen.

    @param data Liste von Messwerten
    @return Arithmetischer Mittelwert
    @throws ValueError Bei leerer Liste
    @lastModified 2026-03-10

    Beispiele:
    >>> mean([1, 2, 3, 4, 5])
    3.0
    >>> mean([10, 20])
    15.0
    """
    if not data:
        raise ValueError("Datenliste darf nicht leer sein")
    # NumPy-vektorisierte Berechnung statt Python-Schleife
    return float(np.mean(data))


def median(data: list[float]) -> float:
    """
    Berechnet den Median (Zentralwert) einer Datenmenge.

    Bei gerader Anzahl: Mittelwert der beiden mittleren Werte.

    @param data Liste von Messwerten
    @return Median
    @throws ValueError Bei leerer Liste
    @lastModified 2026-03-10

    Beispiele:
    >>> median([1, 2, 3])
    2.0
    >>> median([1, 3, 5, 7])
    4.0
    """
    if not data:
        raise ValueError("Datenliste darf nicht leer sein")
    # Sortierte Kopie (Original bleibt unverändert)
    sorted_data = sorted(data)
    n = len(sorted_data)
    mid = n // 2
    if n % 2 == 0:
        return (sorted_data[mid - 1] + sorted_data[mid]) / 2.0
    return float(sorted_data[mid])


def mode(data: list[Any]) -> list[Any]:
    """
    Berechnet den/die Modalwert(e) einer Datenmenge.

    Bei mehreren gleich häufigen Werten werden alle zurückgegeben.

    @param data Liste von Messwerten
    @return Liste der häufigsten Werte
    @throws ValueError Bei leerer Liste
    @lastModified 2026-03-10
    """
    if not data:
        raise ValueError("Datenliste darf nicht leer sein")
    # Häufigkeit jedes Werts zählen
    freq: dict = {}
    for x in data:
        freq[x] = freq.get(x, 0) + 1
    max_freq = max(freq.values())
    # Alle Werte mit maximaler Häufigkeit
    return [k for k, v in freq.items() if v == max_freq]


def variance(data: list[float], ddof: int = 1) -> float:
    """
    Berechnet die Varianz einer Datenmenge.

    Population: ddof=0, Stichprobe: ddof=1 (Bessel-Korrektur).

    Formel: s^2 = (1/(n-ddof)) * sum((x_i - mu)^2)

    Implementierung mit np.var() für vektorisierte Berechnung (kein Python-Loop).

    @param data Liste von Messwerten
    @param ddof Delta Freiheitsgrade (0=Population, 1=Stichprobe)
    @return Varianz
    @throws ValueError Bei zu wenigen Datenpunkten
    @lastModified 2026-03-10

    Beispiele:
    >>> round(variance([2, 4, 4, 4, 5, 5, 7, 9]), 4)
    4.5714
    >>> variance([1, 1, 1], ddof=0)
    0.0
    """
    n = len(data)
    if n <= ddof:
        raise ValueError(f"Zu wenige Datenpunkte (n={n}, ddof={ddof})")
    # NumPy berechnet Varianz direkt vektorisiert: np.var(data, ddof=ddof)
    return float(np.var(data, ddof=ddof))


def std_dev(data: list[float], ddof: int = 1) -> float:
    """
    Berechnet die Standardabweichung einer Datenmenge.

    Formel: s = sqrt(variance(data, ddof))

    Implementierung mit np.std() für vektorisierte Berechnung (kein Python-Loop).

    @param data Liste von Messwerten
    @param ddof Delta Freiheitsgrade (0=Population, 1=Stichprobe)
    @return Standardabweichung
    @lastModified 2026-03-10
    """
    # Validierung über variance() sicherstellen
    n = len(data)
    if n <= ddof:
        raise ValueError(f"Zu wenige Datenpunkte (n={n}, ddof={ddof})")
    # NumPy berechnet Standardabweichung direkt vektorisiert
    return float(np.std(data, ddof=ddof))


def quartiles(data: list[float]) -> tuple[float, float, float]:
    """
    Berechnet die drei Quartile (Q1, Q2, Q3) einer Datenmenge.

    Q2 = Median, Q1 = Median der unteren Hälfte, Q3 = Median der oberen Hälfte.

    @param data Liste von Messwerten
    @return Tupel (Q1, Q2, Q3)
    @throws ValueError Bei zu wenigen Datenpunkten
    @lastModified 2026-03-10
    """
    if len(data) < 4:
        raise ValueError("Mindestens 4 Datenpunkte für Quartile erforderlich")
    sorted_data = sorted(data)
    n = len(sorted_data)

    # Median (Q2)
    q2 = median(sorted_data)

    # Untere und obere Hälfte (ohne Median bei ungerader Anzahl)
    mid = n // 2
    lower = sorted_data[:mid]
    upper = sorted_data[mid + (n % 2):]

    q1 = median(lower)
    q3 = median(upper)
    return q1, q2, q3


def iqr(data: list[float]) -> float:
    """
    Berechnet den Interquartilsabstand (IQR = Q3 - Q1).

    Der IQR ist ein robustes Streuungsmaß (resistent gegen Ausreißer).

    @param data Liste von Messwerten
    @return Interquartilsabstand
    @lastModified 2026-03-10
    """
    q1, _, q3 = quartiles(data)
    return q3 - q1


def skewness(data: list[float]) -> float:
    """
    Berechnet die Schiefe (3. standardisiertes Moment) einer Verteilung.

    Positive Schiefe: rechtsschief (langer rechter Schwanz).
    Negative Schiefe: linksschief (langer linker Schwanz).

    Formel: g1 = (1/n) * sum((x_i - mu)^3) / sigma^3

    Implementierung mit NumPy-Vektorisierung: keine Python-for-Schleife.

    @param data Liste von Messwerten
    @return Fisher-Schiefe
    @throws ValueError Bei weniger als 3 Datenpunkten
    @lastModified 2026-03-10
    """
    n = len(data)
    if n < 3:
        raise ValueError("Mindestens 3 Datenpunkte für Schiefe erforderlich")
    # NumPy-Array für vektorisierte Operationen
    arr = np.asarray(data, dtype=float)
    mu = np.mean(arr)
    sigma = np.std(arr, ddof=0)
    if sigma == 0:
        return 0.0
    # Drittes zentrales Moment vektorisiert: (x - mu)^3 für alle x auf einmal
    m3 = np.mean((arr - mu) ** 3)
    return float(m3 / sigma ** 3)


def kurtosis(data: list[float]) -> float:
    """
    Berechnet die Exzess-Kurtosis (4. Moment - 3) einer Verteilung.

    Normalverteilung: Kurtosis = 0 (Mesokurtisch).
    Positiv (leptokurtisch): spitzere Kurve mit schweren Rändern.
    Negativ (platykurtisch): flachere Kurve.

    Formel: g2 = (1/n) * sum((x_i - mu)^4) / sigma^4 - 3

    Implementierung mit NumPy-Vektorisierung: keine Python-for-Schleife.

    @param data Liste von Messwerten
    @return Exzess-Kurtosis (Fisher-Definition)
    @throws ValueError Bei weniger als 4 Datenpunkten
    @lastModified 2026-03-10
    """
    n = len(data)
    if n < 4:
        raise ValueError("Mindestens 4 Datenpunkte für Kurtosis erforderlich")
    # NumPy-Array für vektorisierte Operationen
    arr = np.asarray(data, dtype=float)
    mu = np.mean(arr)
    sigma = np.std(arr, ddof=0)
    if sigma == 0:
        return 0.0
    # Viertes zentrales Moment vektorisiert: (x - mu)^4 für alle x auf einmal
    m4 = np.mean((arr - mu) ** 4)
    return float(m4 / sigma ** 4 - 3.0)


# ─────────────────────────────────────────────────────────────────────────────
# Normalverteilung
# ─────────────────────────────────────────────────────────────────────────────

def normal_pdf(x: float, mu: float = 0.0, sigma: float = 1.0) -> float:
    """
    Berechnet die Wahrscheinlichkeitsdichtefunktion (PDF) der Normalverteilung.

    Formel: f(x) = (1 / (sigma * sqrt(2*pi))) * exp(-0.5 * ((x-mu)/sigma)^2)

    @param x Stelle, an der die Dichte berechnet wird
    @param mu Erwartungswert (Lageparameter)
    @param sigma Standardabweichung (Skalierungsparameter, sigma > 0)
    @return Wahrscheinlichkeitsdichte f(x)
    @throws ValueError Bei sigma <= 0
    @lastModified 2026-03-10
    """
    if sigma <= 0:
        raise ValueError(f"Standardabweichung sigma muss > 0 sein, nicht {sigma}")
    # Normierungskonstante
    coefficient = 1.0 / (sigma * math.sqrt(2 * math.pi))
    # Exponent des Gauss-Terms
    exponent = -0.5 * ((x - mu) / sigma) ** 2
    return coefficient * math.exp(exponent)


def normal_cdf(x: float, mu: float = 0.0, sigma: float = 1.0) -> float:
    """
    Berechnet die kumulative Verteilungsfunktion (CDF) der Normalverteilung.

    Gibt P(X <= x) für X ~ N(mu, sigma^2) an.

    Formel: F(x) = 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))))

    @param x Obere Integrationsgrenze
    @param mu Erwartungswert
    @param sigma Standardabweichung (sigma > 0)
    @return Wahrscheinlichkeit P(X <= x)
    @lastModified 2026-03-10
    """
    if sigma <= 0:
        raise ValueError(f"sigma muss > 0 sein, nicht {sigma}")
    # Standardisierung und Fehlerfunktion
    z = (x - mu) / (sigma * math.sqrt(2.0))
    return 0.5 * (1.0 + _erf(z))


def normal_ppf(p: float, mu: float = 0.0, sigma: float = 1.0) -> float:
    """
    Berechnet die Quantilfunktion (PPF/Inverse CDF) der Normalverteilung.

    Gibt den Wert x zurück, sodass P(X <= x) = p für X ~ N(mu, sigma^2).

    Formel: x = mu + sigma * sqrt(2) * erfinv(2p - 1)

    @param p Wahrscheinlichkeit (0 < p < 1)
    @param mu Erwartungswert
    @param sigma Standardabweichung (sigma > 0)
    @return x-Wert mit CDF(x) = p
    @throws ValueError Bei p ausserhalb (0,1) oder sigma <= 0
    @lastModified 2026-03-10
    """
    if not 0 < p < 1:
        raise ValueError(f"p muss in (0,1) liegen, nicht {p}")
    if sigma <= 0:
        raise ValueError(f"sigma muss > 0 sein, nicht {sigma}")
    # Inverse der CDF über die inverse Fehlerfunktion
    return mu + sigma * math.sqrt(2.0) * _erfinv(2 * p - 1)


# ─────────────────────────────────────────────────────────────────────────────
# Binomialverteilung
# ─────────────────────────────────────────────────────────────────────────────

def _binom_coeff(n: int, k: int) -> int:
    """
    Berechnet den Binomialkoeffizienten C(n,k) = n! / (k! * (n-k)!).

    @param n Gesamtanzahl
    @param k Auswahlanzahl
    @return Binomialkoeffizient
    @lastModified 2026-03-10
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    # Symmetrie nutzen für Effizienz
    k = min(k, n - k)
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


def binomial_pmf(k: int, n: int, p: float) -> float:
    """
    Berechnet die Wahrscheinlichkeitsmassfunktion (PMF) der Binomialverteilung.

    P(X = k) = C(n,k) * p^k * (1-p)^(n-k)

    @param k Anzahl Erfolge (0 <= k <= n)
    @param n Anzahl Versuche
    @param p Erfolgswahrscheinlichkeit (0 <= p <= 1)
    @return P(X = k)
    @throws ValueError Bei ungültigen Parametern
    @lastModified 2026-03-10
    """
    if n < 0:
        raise ValueError("n muss >= 0 sein")
    if not 0 <= p <= 1:
        raise ValueError("p muss in [0,1] liegen")
    if k < 0 or k > n:
        return 0.0
    # Spezialfälle für p=0 und p=1
    if p == 0:
        return 1.0 if k == 0 else 0.0
    if p == 1:
        return 1.0 if k == n else 0.0
    # Log-Summe für numerische Stabilität bei großen n
    log_prob = (math.log(_binom_coeff(n, k))
                + k * math.log(p)
                + (n - k) * math.log(1 - p))
    return math.exp(log_prob)


def binomial_cdf(k: int, n: int, p: float) -> float:
    """
    Berechnet die kumulative Verteilungsfunktion (CDF) der Binomialverteilung.

    P(X <= k) = sum_{i=0}^{k} C(n,i) * p^i * (1-p)^(n-i)

    @param k Obere Grenze (inklusive)
    @param n Anzahl Versuche
    @param p Erfolgswahrscheinlichkeit (0 <= p <= 1)
    @return P(X <= k)
    @lastModified 2026-03-10
    """
    if k >= n:
        return 1.0
    if k < 0:
        return 0.0
    # Aufsummieren der PMF-Werte
    return sum(binomial_pmf(i, n, p) for i in range(k + 1))


# ─────────────────────────────────────────────────────────────────────────────
# Poisson-Verteilung
# ─────────────────────────────────────────────────────────────────────────────

def poisson_pmf(k: int, lam: float) -> float:
    """
    Berechnet die Wahrscheinlichkeitsmassfunktion (PMF) der Poisson-Verteilung.

    P(X = k) = (lambda^k * e^(-lambda)) / k!

    Die Poisson-Verteilung modelliert die Anzahl seltener Ereignisse
    in einem festen Zeitintervall.

    @param k Anzahl Ereignisse (k >= 0)
    @param lam Erwartungswert lambda (lambda > 0)
    @return P(X = k)
    @throws ValueError Bei ungültigen Parametern
    @lastModified 2026-03-10
    """
    if k < 0:
        return 0.0
    if lam < 0:
        raise ValueError("lambda muss >= 0 sein")
    if lam == 0:
        return 1.0 if k == 0 else 0.0
    # Log-Darstellung für numerische Stabilität
    log_prob = k * math.log(lam) - lam - math.lgamma(k + 1)
    return math.exp(log_prob)


def poisson_cdf(k: int, lam: float) -> float:
    """
    Berechnet die kumulative Verteilungsfunktion (CDF) der Poisson-Verteilung.

    P(X <= k) = sum_{i=0}^{k} (lambda^i * e^(-lambda)) / i!

    @param k Obere Grenze (inklusive)
    @param lam Erwartungswert lambda
    @return P(X <= k)
    @lastModified 2026-03-10
    """
    if k < 0:
        return 0.0
    return sum(poisson_pmf(i, lam) for i in range(k + 1))


# ─────────────────────────────────────────────────────────────────────────────
# Hypothesentests
# ─────────────────────────────────────────────────────────────────────────────

def t_test_one_sample(data: list[float], mu: float) -> tuple[float, float]:
    """
    Führt einen einseitigen Ein-Stichproben-t-Test durch.

    Nullhypothese H0: Populationsmittel = mu
    Alternativhypothese H1: Populationsmittel != mu (zweiseitig)

    Teststatistik: t = (x_bar - mu) / (s / sqrt(n))

    @param data Stichprobendaten
    @param mu Hypothetischer Populationsmittelwert
    @return Tupel (t-Statistik, p-Wert)
    @throws ValueError Bei zu wenigen Datenpunkten
    @lastModified 2026-03-10
    """
    n = len(data)
    if n < 2:
        raise ValueError("Mindestens 2 Datenpunkte erforderlich")
    x_bar = mean(data)
    s = std_dev(data, ddof=1)  # Stichproben-Standardabweichung

    if s == 0:
        # Alle Werte gleich: t-Wert ist 0 wenn mu=x_bar, sonst unendlich
        t_stat = 0.0 if x_bar == mu else float('inf')
        p_val = 1.0 if t_stat == 0.0 else 0.0
        return t_stat, p_val

    # t-Statistik: Abweichung in Einheiten des Standardfehlers
    t_stat = (x_bar - mu) / (s / math.sqrt(n))
    df = n - 1  # Freiheitsgrade

    # Zweiseitiger p-Wert: P(|T| >= |t|) = 2 * P(T >= |t|)
    p_val = 2 * (1 - _t_cdf(abs(t_stat), df))
    return t_stat, p_val


def t_test_two_sample(
    data1: list[float],
    data2: list[float],
    equal_var: bool = True
) -> tuple[float, float]:
    """
    Führt einen Zwei-Stichproben-t-Test (Welch oder Student) durch.

    Nullhypothese H0: mu1 = mu2
    Alternativhypothese H1: mu1 != mu2 (zweiseitig)

    @param data1 Erste Stichprobe
    @param data2 Zweite Stichprobe
    @param equal_var True = Student-t-Test (gleiche Varianzen),
                     False = Welch-t-Test (ungleiche Varianzen)
    @return Tupel (t-Statistik, p-Wert)
    @lastModified 2026-03-10
    """
    n1, n2 = len(data1), len(data2)
    if n1 < 2 or n2 < 2:
        raise ValueError("Jede Stichprobe braucht mindestens 2 Datenpunkte")

    x1, x2 = mean(data1), mean(data2)
    s1, s2 = std_dev(data1, ddof=1), std_dev(data2, ddof=1)

    if equal_var:
        # Gepoolte Standardabweichung (Student)
        sp2 = ((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2)
        se = math.sqrt(sp2 * (1/n1 + 1/n2))
        df = n1 + n2 - 2
    else:
        # Welch-t-Test (keine Annahme gleicher Varianzen)
        se = math.sqrt(s1**2/n1 + s2**2/n2)
        # Welch-Satterthwaite-Freiheitsgrade
        df = int((s1**2/n1 + s2**2/n2)**2 /
                 ((s1**2/n1)**2/(n1-1) + (s2**2/n2)**2/(n2-1)))

    if se == 0:
        t_stat = 0.0 if x1 == x2 else float('inf')
        p_val = 1.0 if t_stat == 0.0 else 0.0
        return t_stat, p_val

    t_stat = (x1 - x2) / se
    p_val = 2 * (1 - _t_cdf(abs(t_stat), df))
    return t_stat, p_val


def chi_square_test(
    observed: list[float],
    expected: list[float]
) -> tuple[float, float]:
    """
    Führt den Chi-Quadrat-Anpassungstest durch.

    Vergleicht beobachtete mit erwarteten Häufigkeiten.
    Nullhypothese H0: Die Daten folgen der erwarteten Verteilung.

    Teststatistik: chi2 = sum((O_i - E_i)^2 / E_i)

    @param observed Beobachtete Häufigkeiten
    @param expected Erwartete Häufigkeiten
    @return Tupel (chi2-Statistik, p-Wert)
    @throws ValueError Bei ungleicher Länge oder negativen erwarteten Werten
    @lastModified 2026-03-10
    """
    if len(observed) != len(expected):
        raise ValueError("observed und expected müssen gleich lang sein")
    if any(e <= 0 for e in expected):
        raise ValueError("Alle erwarteten Häufigkeiten müssen > 0 sein")

    # Chi-Quadrat-Statistik: quadratische Abweichungen normiert durch Erwartung
    chi2 = sum((o - e)**2 / e for o, e in zip(observed, expected))
    # Freiheitsgrade: Anzahl Kategorien - 1
    df = len(observed) - 1

    # p-Wert: P(Chi2 >= chi2_stat)
    p_val = 1.0 - _chi2_cdf(chi2, df)
    return chi2, p_val


# ─────────────────────────────────────────────────────────────────────────────
# Bayes-Theorem
# ─────────────────────────────────────────────────────────────────────────────

def bayes_theorem(prior: float, likelihood: float, evidence: float) -> float:
    """
    Berechnet die A-posteriori-Wahrscheinlichkeit nach dem Bayes-Theorem.

    Formel: P(A|B) = P(B|A) * P(A) / P(B)

    Dabei sind:
    - prior     = P(A)        : A-priori-Wahrscheinlichkeit der Hypothese
    - likelihood = P(B|A)     : Likelihood der Beobachtung gegeben die Hypothese
    - evidence  = P(B)        : Gesamtwahrscheinlichkeit der Beobachtung

    @param prior A-priori-Wahrscheinlichkeit P(A)
    @param likelihood Bedingte Wahrscheinlichkeit P(B|A)
    @param evidence Marginale Wahrscheinlichkeit P(B)
    @return A-posteriori-Wahrscheinlichkeit P(A|B)
    @throws ValueError Bei evidence = 0
    @lastModified 2026-03-10
    """
    if evidence == 0:
        raise ValueError("Evidence P(B) darf nicht 0 sein")
    return likelihood * prior / evidence


# ─────────────────────────────────────────────────────────────────────────────
# Monte-Carlo-Simulation
# ─────────────────────────────────────────────────────────────────────────────

def monte_carlo_pi(n_samples: int = 100000, seed: int = None) -> float:
    """
    Schätzt die Kreiszahl pi durch Monte-Carlo-Simulation.

    Methode: Zufällige Punkte im Einheitsquadrat [0,1]^2 werden generiert.
    Der Anteil der Punkte innerhalb des Einheitskreisviertels (x^2 + y^2 <= 1)
    approximiert pi/4.

    Formel: pi ≈ 4 * (Punkte im Kreis) / (Gesamtpunkte)

    Konvergenzrate: O(1/sqrt(n)) durch das Gesetz der großen Zahlen.

    @param n_samples Anzahl zufälliger Stichproben
    @param seed Zufallsseed für Reproduzierbarkeit (None = zufällig)
    @return Näherung von pi
    @throws ValueError Bei n_samples <= 0
    @lastModified 2026-03-10
    """
    if n_samples <= 0:
        raise ValueError("n_samples muss > 0 sein")

    # Zufallsgenerator initialisieren
    rng = random.Random(seed)

    # Zähler für Punkte innerhalb des Viertelkreises
    inside = 0
    for _ in range(n_samples):
        x = rng.random()
        y = rng.random()
        # Punkt liegt im Viertelkreis wenn Abstand vom Ursprung <= 1
        if x * x + y * y <= 1.0:
            inside += 1

    # pi ≈ 4 * (Anteil der Punkte im Kreis)
    return 4.0 * inside / n_samples
