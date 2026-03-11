"""
@file arbitrary_precision.py
@brief Beliebig genaue Arithmetik mit mpmath: ζ-Funktion, π, e, Gamma, Bernoulli-Zahlen.
@description
    Dieses Modul stellt Funktionen für beliebig genaue numerische Berechnungen
    bereit, die mpmath (Multiple Precision Math) als Backend nutzen.

    Funktionen:
    - set_precision()           – Setzt globale Rechenpräzision
    - riemann_zeta_highprec()   – ζ(s) mit beliebiger Genauigkeit
    - riemann_zero_verify()     – Verifikation der n-ten Riemann-Nullstelle
    - pi_highprec()             – π auf beliebig viele Dezimalstellen
    - e_highprec()              – e (Euler-Zahl) auf beliebig viele Dezimalstellen
    - gamma_highprec()          – Γ(s) mit hoher Präzision
    - bernoulli_highprec()      – Bernoulli-Zahl B_n mit hoher Präzision
    - continued_fraction_expansion() – Kettenbruchentwicklung einer Zahl

    Mathematische Grundlagen:
    - mpmath arbeitet mit beliebiger Dezimalstellenanzahl (mp.dps)
    - Riemann-Hypothese: Alle nicht-trivialen Nullstellen von zeta(s) liegen auf Re(s) = 1/2
    - Bernoulli-Zahlen: B_n = (-1)^n * n! * Summe_{k=0}^{n} (-1)^k * C(n,k) / (k+1)
    - Kettenbruchentwicklung: x = a_0 + 1/(a_1 + 1/(a_2 + ...))

@author Michael Fuhrmann
@lastModified 2026-03-10
"""

import mpmath
from typing import Union


# =============================================================================
# PRÄZISIONS-MANAGEMENT
# =============================================================================

def set_precision(decimal_digits: int) -> None:
    """
    @brief Setzt die globale Rechenpräzision in Dezimalstellen.
    @description
        mpmath.mp.dps (decimal places) steuert die Anzahl der signifikanten
        Dezimalstellen für alle mpmath-Berechnungen.

        Standardwert: 15 (entspricht float64-Genauigkeit)
        Empfehlung für Riemann-Nullstellen: ≥ 50
        Empfehlung für π/e: ≥ 1000 für sichtbaren Effekt
    @param decimal_digits Anzahl der Dezimalstellen (muss > 0 sein)
    @lastModified 2026-03-10
    """
    # Setze globale Präzision in mpmath
    mpmath.mp.dps = decimal_digits


# =============================================================================
# RIEMANN-ZETA-FUNKTION
# =============================================================================

def riemann_zeta_highprec(s: complex, digits: int = 50) -> complex:
    """
    @brief Berechnet ζ(s) mit beliebiger Genauigkeit via mpmath.
    @description
        Die Riemann-Zeta-Funktion ist definiert als:
        ζ(s) = Σ_{n=1}^{∞} n^{-s}  für Re(s) > 1

        Durch analytische Fortsetzung gilt sie für alle s ∈ ℂ außer s = 1
        (Pol erster Ordnung mit Residuum 1).

        Spezielle Werte:
        - ζ(2) = π²/6 ≈ 1.6449...
        - ζ(4) = π⁴/90
        - ζ(-1) = -1/12  (regularisiert)

        mpmath nutzt die Euler-Maclaurin-Formel für Re(s) > 0 und
        die Funktionalgleichung für Re(s) ≤ 0.
    @param s      Komplexes Argument s aus C ohne {1}
    @param digits Anzahl der Dezimalstellen (Standard: 50)
    @return ζ(s) als komplexe Zahl
    @lastModified 2026-03-10
    """
    # Präzision setzen
    mpmath.mp.dps = digits
    # Berechne ζ(s) mit mpmath
    result = mpmath.zeta(s)
    return complex(result)


def riemann_zero_verify(n: int, digits: int = 100) -> dict:
    """
    @brief Verifiziert die n-te nicht-triviale Riemann-Nullstelle.
    @description
        Die Riemann-Hypothese besagt, dass alle nicht-trivialen Nullstellen
        ρ_n = 1/2 + i·t_n von ζ(s) auf der kritischen Linie Re(s) = 1/2 liegen.

        mpmath.zetazero(n) berechnet die n-te Nullstelle numerisch auf
        beliebig viele Stellen.

        Bekannte Werte:
        - ρ₁ ≈ 1/2 + 14.134725... · i
        - ρ₂ ≈ 1/2 + 21.022040... · i
        - ρ₃ ≈ 1/2 + 25.010858... · i

        Bisher wurden alle bekannten Nullstellen auf der kritischen Linie
        numerisch verifiziert.
    @param n      Index der Nullstelle (n ≥ 1)
    @param digits Anzahl der Präzisionsstellen (Standard: 100)
    @return Dictionary mit Informationen über die Nullstelle
    @raises ValueError falls n < 1
    @lastModified 2026-03-10
    """
    # Validiere Index
    if n < 1:
        raise ValueError(f"Index n muss ≥ 1 sein, erhalten: {n}")

    # Präzision setzen
    mpmath.mp.dps = digits

    # Berechne n-te Nullstelle
    zero = mpmath.zetazero(n)

    # Extrahiere Real- und Imaginärteil
    re_part = float(mpmath.re(zero))
    im_part = mpmath.im(zero)

    # Prüfe ob Nullstelle auf der kritischen Linie Re(ρ) = 1/2 liegt
    # Toleranz: 10^(-digits/2)
    tolerance = mpmath.mpf(10) ** (-(digits // 2))
    on_critical_line = abs(mpmath.re(zero) - mpmath.mpf('0.5')) < tolerance

    return {
        'n': n,
        'real_part': re_part,
        'imag_part': mpmath.nstr(im_part, digits),
        'on_critical_line': on_critical_line,
        'zero_value': complex(zero)
    }


# =============================================================================
# KONSTANTEN MIT HOHER PRÄZISION
# =============================================================================

def pi_highprec(digits: int = 1000) -> str:
    """
    @brief Berechnet π auf die gewünschte Anzahl Dezimalstellen.
    @description
        Kreiszahl π = 3.14159265358979323846...

        mpmath nutzt den Chudnovsky-Algorithmus für schnelle π-Berechnung:
        1/π = (12 / 640320^{3/2}) · Σ (-1)^k (6k)! (13591409 + 545140134k) /
              ((3k)! (k!)^3 640320^{3k})

        Diese Formel konvergiert mit ~14 Dezimalstellen pro Term.
    @param digits Anzahl der Dezimalstellen (Standard: 1000)
    @return π als Dezimalstring mit 'digits' Stellen
    @lastModified 2026-03-10
    """
    # Präzision setzen (etwas mehr für Rundungssicherheit)
    mpmath.mp.dps = digits + 10

    # Berechne π mit mpmath
    pi_value = mpmath.pi

    # Ausgabe als String mit gewünschter Stellenanzahl
    return mpmath.nstr(pi_value, digits)


def e_highprec(digits: int = 1000) -> str:
    """
    @brief Berechnet die Euler-Zahl e auf die gewünschte Anzahl Dezimalstellen.
    @description
        Eulersche Zahl e = 2.71828182845904523536...

        Berechnung via:
        e = Σ_{n=0}^{∞} 1/n! = 1 + 1 + 1/2 + 1/6 + 1/24 + ...

        mpmath nutzt eine effiziente binäre Aufteilungsformel.
    @param digits Anzahl der Dezimalstellen (Standard: 1000)
    @return e als Dezimalstring mit 'digits' Stellen
    @lastModified 2026-03-10
    """
    # Präzision setzen
    mpmath.mp.dps = digits + 10

    # Berechne e mit mpmath
    e_value = mpmath.e

    # Ausgabe als String mit gewünschter Stellenanzahl
    return mpmath.nstr(e_value, digits)


# =============================================================================
# GAMMA-FUNKTION
# =============================================================================

def gamma_highprec(s: complex, digits: int = 50) -> complex:
    """
    @brief Berechnet die Gamma-Funktion Γ(s) mit hoher Präzision.
    @description
        Die Gamma-Funktion verallgemeinert die Fakultät:
        Γ(n) = (n-1)!  für positive ganze Zahlen n
        Γ(s) = ∫₀^∞ t^{s-1} e^{-t} dt  für Re(s) > 0

        Durch die Rekursionsformel Gamma(s+1) = s*Gamma(s) kann sie auf C ohne {0,-1,-2,...}
        analytisch fortgesetzt werden.

        Besondere Werte:
        - Γ(1/2) = √π
        - Γ(1) = 1
        - Γ(n) = (n-1)! für n ∈ ℕ
    @param s      Komplexes Argument (nicht 0, -1, -2, ...)
    @param digits Anzahl der Dezimalstellen (Standard: 50)
    @return Γ(s) als komplexe Zahl
    @lastModified 2026-03-10
    """
    # Präzision setzen
    mpmath.mp.dps = digits

    # Berechne Γ(s) mit mpmath
    result = mpmath.gamma(s)
    return complex(result)


# =============================================================================
# BERNOULLI-ZAHLEN
# =============================================================================

def bernoulli_highprec(n: int, digits: int = 50) -> str:
    """
    @brief Berechnet die n-te Bernoulli-Zahl B_n mit hoher Präzision.
    @description
        Die Bernoulli-Zahlen sind rationale Zahlen, die über die erzeugende
        Funktion definiert sind:
        t / (e^t - 1) = Σ_{n=0}^{∞} B_n · t^n / n!

        Erste Werte:
        - B_0 = 1
        - B_1 = -1/2
        - B_2 = 1/6
        - B_3 = 0
        - B_4 = -1/30

        Für ungerade n > 1 gilt: B_n = 0

        Zusammenhang mit der Riemann-Zeta-Funktion:
        ζ(-n) = -B_{n+1} / (n+1)  für n ≥ 0
    @param n      Index der Bernoulli-Zahl (n ≥ 0)
    @param digits Anzahl der Dezimalstellen (Standard: 50)
    @return B_n als Dezimalstring
    @raises ValueError falls n < 0
    @lastModified 2026-03-10
    """
    # Validiere Index
    if n < 0:
        raise ValueError(f"Index n muss ≥ 0 sein, erhalten: {n}")

    # Präzision setzen
    mpmath.mp.dps = digits + 10

    # Berechne Bernoulli-Zahl mit mpmath
    bn = mpmath.bernoulli(n)

    # Ausgabe als Dezimalstring
    return mpmath.nstr(bn, digits)


# =============================================================================
# KETTENBRUCHENTWICKLUNG
# =============================================================================

def continued_fraction_expansion(x: float, n_terms: int = 20) -> list:
    """
    @brief Berechnet die Kettenbruchentwicklung einer reellen Zahl.
    @description
        Jede reelle Zahl x lässt sich als Kettenbruch darstellen:
        x = a_0 + 1/(a_1 + 1/(a_2 + 1/(a_3 + ...)))

        Kurzschreibweise: x = [a_0; a_1, a_2, a_3, ...]

        wobei a_i = floor(x_i) und x_{i+1} = 1 / (x_i - a_i)

        Beispiele:
        - π = [3; 7, 15, 1, 292, 1, 1, 1, 2, ...]
        - e = [2; 1, 2, 1, 1, 4, 1, 1, 6, ...]
        - √2 = [1; 2, 2, 2, 2, ...] (periodisch)
        - Goldener Schnitt φ = [1; 1, 1, 1, ...] (alle Einsen)

        Rationale Zahlen haben endliche Kettenbrüche.
        Irrationale Zahlen haben unendliche Kettenbrüche.
        Quadratische Irrationalitäten haben periodische Kettenbrüche.
    @param x       Die zu entwickelnde reelle Zahl
    @param n_terms Maximale Anzahl der Terme (Standard: 20)
    @return Liste der Kettenbruchkoeffizienten [a_0, a_1, a_2, ...]
    @lastModified 2026-03-10
    """
    # Verwende mpmath für genaue Berechnung
    mpmath.mp.dps = 50

    # Koeffizientenliste initialisieren
    coefficients = []

    # Konvertiere zu mpmath-Zahl für Präzision
    x_mp = mpmath.mpf(str(x))

    for _ in range(n_terms):
        # a_i = floor(x_i) (ganzzahliger Anteil)
        a_i = int(mpmath.floor(x_mp))
        coefficients.append(a_i)

        # Berechne Rest: x_{i+1} = 1 / (x_i - a_i)
        remainder = x_mp - a_i

        # Beende wenn Rest sehr klein (rationale Zahl oder Präzisionsgrenze)
        if abs(remainder) < mpmath.mpf(10) ** (-40):
            break

        # Nächster Wert
        x_mp = mpmath.mpf(1) / remainder

    return coefficients


# =============================================================================
# ERWEITERTE HILFSFUNKTIONEN
# =============================================================================

def euler_mascheroni_highprec(digits: int = 50) -> str:
    """
    @brief Berechnet die Euler-Mascheroni-Konstante γ mit hoher Präzision.
    @description
        Die Euler-Mascheroni-Konstante γ ist definiert als:
        γ = lim_{n→∞} (Σ_{k=1}^{n} 1/k - ln(n)) ≈ 0.5772156649...

        Sie tritt in vielen Formeln der Analysis und Zahlentheorie auf,
        z.B. im Zusammenhang mit der Gamma-Funktion:
        Γ'(1) = -γ

        Es ist unbekannt, ob γ irrational ist.
    @param digits Anzahl der Dezimalstellen (Standard: 50)
    @return γ als Dezimalstring
    @lastModified 2026-03-10
    """
    # Präzision setzen
    mpmath.mp.dps = digits + 10

    # Berechne Euler-Mascheroni-Konstante
    gamma_const = mpmath.euler

    return mpmath.nstr(gamma_const, digits)


def log_highprec(x: float, base: float = None, digits: int = 50) -> str:
    """
    @brief Berechnet den Logarithmus mit hoher Präzision.
    @description
        Berechnet ln(x) (natürlicher Logarithmus) oder log_b(x) (Logarithmus
        zur Basis b) mit beliebiger Genauigkeit.

        Mathematisch:
        - ln(x) = ∫₁^x dt/t
        - log_b(x) = ln(x) / ln(b)
    @param x      Das Argument (x > 0)
    @param base   Basis des Logarithmus (None = natürlicher Log)
    @param digits Anzahl der Dezimalstellen (Standard: 50)
    @return Logarithmus als Dezimalstring
    @raises ValueError falls x ≤ 0
    @lastModified 2026-03-10
    """
    # Validierung
    if x <= 0:
        raise ValueError(f"Argument x muss positiv sein, erhalten: {x}")

    # Präzision setzen
    mpmath.mp.dps = digits + 10

    # Berechne Logarithmus
    if base is None:
        result = mpmath.log(x)
    else:
        result = mpmath.log(x, base)

    return mpmath.nstr(result, digits)


def sqrt_highprec(x: float, digits: int = 50) -> str:
    """
    @brief Berechnet √x mit hoher Präzision.
    @description
        Berechnet die Quadratwurzel einer nicht-negativen reellen Zahl
        mit beliebiger Genauigkeit.

        mpmath nutzt Newton-Iteration mit automatischer Schrittweitensteuerung.
    @param x      Das Argument (x ≥ 0)
    @param digits Anzahl der Dezimalstellen (Standard: 50)
    @return √x als Dezimalstring
    @raises ValueError falls x < 0
    @lastModified 2026-03-10
    """
    # Validierung
    if x < 0:
        raise ValueError(f"Argument x muss nicht-negativ sein, erhalten: {x}")

    # Präzision setzen
    mpmath.mp.dps = digits + 10

    # Berechne Quadratwurzel
    result = mpmath.sqrt(x)

    return mpmath.nstr(result, digits)


# =============================================================================
# RIEMANN-NULLSTELLEN (mpmath)
# =============================================================================

def zeta_zeros_mpmath(n: int, prec: int = 50) -> list:
    """
    @brief Berechnet die ersten n Riemann-Nullstellen mit prec Dezimalstellen.
    @description
        Die Riemann-Hypothese besagt, dass alle nicht-trivialen Nullstellen
        der Riemann-Zeta-Funktion ζ(s) auf der kritischen Geraden Re(s) = 1/2
        liegen:
            s_k = 1/2 + i·t_k,  t_k ∈ ℝ, t_k > 0

        Die ersten bekannten Nullstellen t_k:
            t_1 ≈ 14.1347...
            t_2 ≈ 21.0220...
            t_3 ≈ 25.0109...
            ...

        mpmath.zetazero(k) berechnet die k-te Nullstelle mit beliebiger
        Genauigkeit.

    @param n    Anzahl der zu berechnenden Nullstellen (muss ≥ 1 sein)
    @param prec Anzahl der Dezimalstellen (Standard: 50)
    @return Liste von mpmath-Zahlen: die ersten n Nullstellen 1/2 + i·t_k
    @raises ValueError falls n < 1 oder prec < 1
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Eingabevalidierung
    if n < 1:
        raise ValueError(f"n muss mindestens 1 sein, erhalten: {n}")
    if prec < 1:
        raise ValueError(f"prec muss mindestens 1 sein, erhalten: {prec}")

    # Präzision für alle Berechnungen setzen
    mpmath.mp.dps = prec + 10  # Puffer für numerische Fehler

    zeros = []
    for k in range(1, n + 1):
        # mpmath.zetazero(k): k-te Nullstelle der Riemann-Zeta-Funktion
        # Gibt komplexe Zahl 1/2 + i·t_k zurück
        zero = mpmath.zetazero(k)
        zeros.append(zero)

    return zeros


def verify_riemann_hypothesis_mpmath(n_zeros: int, prec: int = 50) -> dict:
    """
    @brief Prüft Re(s) = 1/2 für die ersten n Nullstellen der Riemann-Zeta-Funktion.
    @description
        Die Riemann-Hypothese besagt: Alle nicht-trivialen Nullstellen von ζ(s)
        haben Realteil 1/2. Diese Funktion prüft dies numerisch für die
        ersten n_zeros Nullstellen mit prec Dezimalstellen Genauigkeit.

        Vorgehen:
        1. Berechne die k-te Nullstelle s_k = σ_k + i·t_k via zetazero(k)
        2. Prüfe ob |σ_k - 1/2| < 10^{-prec/2}
        3. Sammle Statistiken über alle Nullstellen

        Ergebnis: Numerische Verifikation (kein mathematischer Beweis!).
        Die Riemann-Hypothese ist noch unbewiesen.

    @param n_zeros Anzahl der zu prüfenden Nullstellen (Standard: min. 1)
    @param prec    Genauigkeit in Dezimalstellen (Standard: 50)
    @return Dictionary mit Verifikationsergebnissen:
            - 'n_checked': Anzahl geprüfter Nullstellen
            - 'all_on_critical_line': True wenn alle auf Re=1/2
            - 'max_deviation': Maximale Abweichung von 1/2
            - 'zeros': Liste der Nullstellen als Strings
            - 'real_parts': Realteile der Nullstellen
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Eingabevalidierung
    if n_zeros < 1:
        raise ValueError(f"n_zeros muss mindestens 1 sein, erhalten: {n_zeros}")

    # Präzision setzen (Puffer für Rundungsfehler)
    mpmath.mp.dps = prec + 15

    # Toleranz: Abweichung kleiner als 10^{-(prec/2)} gilt als "auf der Geraden"
    tolerance = mpmath.mpf(10) ** (-(prec // 2))

    # Nullstellen berechnen und prüfen
    zeros_list = zeta_zeros_mpmath(n_zeros, prec=prec + 10)

    real_parts = []
    deviations = []
    max_deviation = mpmath.mpf(0)

    for s_k in zeros_list:
        # Realteil der Nullstelle
        sigma_k = mpmath.re(s_k)
        real_parts.append(float(sigma_k))

        # Abweichung von 1/2
        deviation = abs(sigma_k - mpmath.mpf('0.5'))
        deviations.append(float(deviation))

        # Maximum aktualisieren
        if deviation > max_deviation:
            max_deviation = deviation

    # Prüfung: Alle Realteile innerhalb der Toleranz?
    all_on_line = all(d < float(tolerance) for d in deviations)

    return {
        'n_checked': n_zeros,
        'all_on_critical_line': all_on_line,
        'max_deviation': float(max_deviation),
        'tolerance': float(tolerance),
        'zeros': [mpmath.nstr(z, min(prec, 20)) for z in zeros_list],
        'real_parts': real_parts,
        'imaginary_parts': [float(mpmath.im(z)) for z in zeros_list],
        'hypothesis_consistent': all_on_line  # Numerisch konsistent (kein Beweis!)
    }


def pi_mpmath(prec: int = 100) -> str:
    """
    @brief Berechnet π auf prec Dezimalstellen via mpmath.
    @description
        mpmath verwendet den Chudnovsky-Algorithmus für sehr hohe Präzisionen
        sowie AGM-basierte Methoden (Arithmetic-Geometric Mean).

        Der Chudnovsky-Algorithmus hat Komplexität O(n·(log n)^3) für
        n Dezimalstellen und ist damit einer der schnellsten bekannten.

        Mathematische Identität hinter Chudnovsky:
            1/π = (12/√(640320³)) · Σ_{k=0}^∞ (6k)!·(13591409+545140134k)
                                                 / ((3k)!(k!)^3·(-640320³)^k)

        Bekannte Werte:
            π ≈ 3.14159265358979323846...
            π auf 100 Stellen: 3.14159265358979323846...28841971693993751...

    @param prec Anzahl der Dezimalstellen (Standard: 100)
    @return π als Dezimalstring mit prec Stellen
    @raises ValueError falls prec < 1
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Eingabevalidierung
    if prec < 1:
        raise ValueError(f"prec muss mindestens 1 sein, erhalten: {prec}")

    # Präzision setzen (Puffer für Rundungsfehler beim Konvertieren)
    mpmath.mp.dps = prec + 10

    # π berechnen (mpmath nutzt intern den besten verfügbaren Algorithmus)
    pi_val = mpmath.pi

    # Als Dezimalstring mit genau prec signifikanten Stellen ausgeben
    return mpmath.nstr(pi_val, prec)
