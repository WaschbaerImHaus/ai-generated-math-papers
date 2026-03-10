"""
@file __init__.py
@brief Paket-Initialisierung für das specialist-maths Projekt.
@description
    Diese Datei macht das Verzeichnis src/ zu einem Python-Package und
    stellt alle wichtigen Klassen und Funktionen aus den Submodulen
    direkt unter dem Paketnamen bereit.

    Importierte Module:
    - algebra:             Polynome, Gleichungslöser, Zahlentheorie, RSA
    - analysis:            Ableitung, Integration, Nullstellensuche, Taylor
    - linear_algebra:      Vektoren, Matrizen, Zerlegungen
    - statistics_math:     Deskriptive Statistik, Verteilungen
    - ode:                 Gewöhnliche Differentialgleichungen
    - fourier:             Fourier-Transformation (DFT/FFT/IFFT)
    - complex_analysis:    Gamma-Funktion, Riemann-Zeta, Xi-Funktion
    - numerical_methods:   Interpolation, Optimierung
    - modular_forms:       Modulare Gruppe, Eisenstein-Reihen, J-Invariante
    - p_adic:              P-adische Zahlen und Normen
    - visualization:       2D-Funktionsgraphen
    - proof_theory:        Collatz, Goldbach, Primzahlsätze

    Verwendungsbeispiel:
        from src import Polynomial, Vector, Matrix, riemann_zeta

    Hinweis: safe_parse_expr ist intern als _safe_parse in analysis.py
    implementiert und wird hier als Alias exportiert.

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

# ---------------------------------------------------------------------------
# Algebra-Modul: Polynome, Gleichungslöser, Zahlentheorie, Kryptographie
# ---------------------------------------------------------------------------
from .algebra import (
    Polynomial,           # Polynomklasse mit add/mul/evaluate/differentiate
    solve_linear,         # Löst lineare Gleichung ax + b = 0
    solve_quadratic,      # Löst quadratische Gleichung ax² + bx + c = 0
    gcd,                  # Größter gemeinsamer Teiler (Euklidischer Algorithmus)
    lcm,                  # Kleinstes gemeinsames Vielfaches
    extended_gcd,         # Erweiterter euklidischer Algorithmus (Bezout-Koeffizienten)
    mod_inverse,          # Modulare Inverse a^{-1} mod m
    is_prime,             # Primzahltest (Miller-Rabin, deterministisch)
    prime_factorization,  # Primfaktorzerlegung
    euler_phi,            # Euler'sche Phi-Funktion φ(n)
    rsa_keygen,           # RSA-Schlüsselgenerierung (p, q → public/private key)
    rsa_encrypt,          # RSA-Verschlüsselung: m^e mod n
    rsa_decrypt,          # RSA-Entschlüsselung: c^d mod n
)

# ---------------------------------------------------------------------------
# Analysis-Modul: Differentiation, Integration, Nullstellensuche, Taylor
# ---------------------------------------------------------------------------
from .analysis import (
    _safe_parse as safe_parse_expr,  # Sicheres Parsen: kein eval(), Whitelist-basiert
    numerical_derivative,            # Zentrale Differenzenquotienten (1. und 2. Ordnung)
    numerical_integral,              # Simpson-Regel für bestimmte Integrale
    newton_raphson,                  # Newton-Raphson-Verfahren zur Nullstellensuche
    bisection,                       # Bisektionsverfahren (garantiert konvergent)
    taylor_series,                   # Symbolische Taylor-Reihenentwicklung via SymPy
)

# ---------------------------------------------------------------------------
# Lineare Algebra: Vektoren, Matrizen, numerische Zerlegungen
# ---------------------------------------------------------------------------
from .linear_algebra import (
    Vector,            # Vektorklasse (dot/cross/norm/normalize)
    Matrix,            # Matrixklasse (det/inv/solve/eigenvalues/eigenvectors)
    gram_schmidt,      # Gram-Schmidt-Orthogonalisierung
    lu_decomposition,  # LU-Zerlegung nach Doolittle mit Teilpivotisierung
    qr_decomposition,  # QR-Zerlegung via Householder-Reflexionen
    svd,               # Singulärwertzerlegung A = U Σ Vᵀ
)

# ---------------------------------------------------------------------------
# Statistik-Modul: Deskriptive Statistik, Wahrscheinlichkeitsverteilungen
# ---------------------------------------------------------------------------
from .statistics_math import (
    mean,      # Arithmetisches Mittel
    median,    # Median (50. Perzentil)
    mode,      # Modus (häufigster Wert)
    variance,  # Varianz (mit Bessel-Korrektur n-1)
    std_dev,   # Standardabweichung
)

# ---------------------------------------------------------------------------
# ODE-Modul: Gewöhnliche Differentialgleichungen
# ---------------------------------------------------------------------------
from .ode import (
    euler_method,  # Euler-Vorwärtsverfahren (einfachstes Einschrittverfahren)
    runge_kutta4,  # Klassisches Runge-Kutta-Verfahren 4. Ordnung
)

# ---------------------------------------------------------------------------
# Fourier-Modul: Diskrete und schnelle Fourier-Transformation
# ---------------------------------------------------------------------------
from .fourier import (
    dft,   # Diskrete Fourier-Transformation O(n²)
    fft,   # Schnelle Fourier-Transformation Cooley-Tukey O(n log n)
    ifft,  # Inverse FFT
)

# ---------------------------------------------------------------------------
# Komplexe Analysis: Gamma-Funktion, Riemann-Zeta, Xi-Funktion
# ---------------------------------------------------------------------------
from .complex_analysis import (
    gamma_lanczos,  # Gamma-Funktion via Lanczos-Approximation
    riemann_zeta,   # Riemann-Zeta-Funktion ζ(s) für alle s ∈ ℂ
    xi_function,    # Xi-Funktion ξ(s) = (1/2)s(s-1)π^{-s/2}Γ(s/2)ζ(s)
)

# ---------------------------------------------------------------------------
# Numerische Methoden: Interpolation, Optimierung
# ---------------------------------------------------------------------------
from .numerical_methods import (
    lagrange_interpolation,  # Lagrange-Interpolationspolynom (Funktion)
    NewtonInterpolation,     # Newton-Interpolation mit dividierten Differenzen (Klasse)
)

# ---------------------------------------------------------------------------
# Modulformen: Modulare Gruppe, Eisenstein-Reihen, J-Invariante
# ---------------------------------------------------------------------------
from .modular_forms import (
    ModularGroup,       # SL(2,ℤ) - modulare Gruppe mit Möbius-Transformationen
    eisenstein_series,  # Eisenstein-Reihen E_k(τ) für gerades k ≥ 4
    delta_function,     # Ramanujan-Delta-Funktion Δ(τ) = η(τ)^24
    j_invariant,        # J-Invariante j(τ) = 1728 · E₄³ / (E₄³ - E₆²)
)

# ---------------------------------------------------------------------------
# P-adische Zahlen: P-adische Bewertung, Norm und Arithmetik
# ---------------------------------------------------------------------------
from .p_adic import (
    PAdicNumber,      # Klasse für p-adische Zahlen mit Arithmetik
    p_adic_valuation,  # P-adische Bewertung v_p(n)
    p_adic_norm,      # P-adische Norm |n|_p = p^{-v_p(n)}
)

# ---------------------------------------------------------------------------
# Visualisierung: 2D-Funktionsgraphen
# ---------------------------------------------------------------------------
from .visualization import (
    plot_function_2d,  # 2D-Graph einer Funktion im angegebenen Intervall
)

# ---------------------------------------------------------------------------
# Beweistheorie: Zahlentheoretische Folgen und Vermutungen
# ---------------------------------------------------------------------------
from .proof_theory import (
    collatz_sequence,        # Collatz-Folge ab Startwert n
    goldbach_decomposition,  # Goldbach-Zerlegung: gerade Zahl als Summe zweier Primzahlen
)

# ---------------------------------------------------------------------------
# Versionsinformation
# ---------------------------------------------------------------------------
__version__: str = "11.0"
__author__: str = "Kurt Ingwer"
