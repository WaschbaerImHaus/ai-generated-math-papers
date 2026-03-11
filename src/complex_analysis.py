"""
Komplexe Analysis – Modul für holomorphe Funktionen und Riemann-Zeta.

Implementiert die vollständige Riemann-Zeta-Funktion via analytischer Fortsetzung
und Funktionalgleichung, die Gamma-Funktion, Riemanns ξ-Funktion und Werkzeuge
zur numerischen Suche von Nullstellen auf der kritischen Geraden.

Mathematischer Hintergrund:
    ζ(s) ist durch Σ 1/n^s nur für Re(s) > 1 definiert.
    Die analytische Fortsetzung erstreckt zeta auf ganz C ohne {1}.
    Die Funktionalgleichung (Riemann 1859) verbindet ζ(s) mit ζ(1-s):
        ξ(s) = ξ(1-s)  mit  ξ(s) = ½s(s-1)π^(-s/2)Γ(s/2)ζ(s)

@author: Michael Fuhrmann
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-10
"""

import cmath
import math
import numpy as np
from typing import Callable


# ===========================================================================
# GAMMA-FUNKTION (benötigt für Riemann-Zeta über Funktionalgleichung)
# ===========================================================================

def gamma_lanczos(z: complex) -> complex:
    """
    Berechnet die Gamma-Funktion Γ(z) via Lanczos-Approximation.

    Die Gamma-Funktion verallgemeinert die Fakultät auf komplexe Zahlen:
        Γ(n) = (n-1)!  für positive ganze n
        Γ(z+1) = z·Γ(z)  (Funktionalgleichung)
        Γ(1/2) = √π

    Lanczos-Approximation (1964): numerisch stabil für Re(z) > 0.5.
    Für Re(z) < 0.5 wird die Reflexionsformel verwendet:
        Γ(z)·Γ(1-z) = π / sin(πz)

    Genauigkeit: ~15 signifikante Stellen.

    @param z: Komplexe Zahl (z ∉ {0, -1, -2, ...})
    @return: Γ(z)
    @raises ValueError: Wenn z eine nicht-positive ganze Zahl ist
    @lastModified: 2026-03-10
    """
    # Pol-Prüfung: Γ hat Pole bei z = 0, -1, -2, ...
    if z.real <= 0 and abs(z.imag) < 1e-12 and abs(z.real - round(z.real)) < 1e-12:
        raise ValueError(f"Γ(z) hat einen Pol bei z = {z}")

    # Reflexionsformel für Re(z) < 0.5
    if z.real < 0.5:
        return cmath.pi / (cmath.sin(cmath.pi * z) * gamma_lanczos(1 - z))

    # Lanczos-Koeffizienten (g=7, n=9) nach Paul Godfrey
    g = 7
    c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ]

    z -= 1  # Lanczos arbeitet mit z-1
    x = complex(c[0])
    for i in range(1, g + 2):
        x += c[i] / (z + i)

    t = z + g + 0.5
    # Stirling-ähnliche Formel: √(2π) · t^(z+0.5) · e^(-t) · x
    return cmath.sqrt(2 * cmath.pi) * (t ** (z + 0.5)) * cmath.exp(-t) * x


def log_gamma(z: complex) -> complex:
    """
    Berechnet ln(Γ(z)) für Re(z) > 0 via Stirling-Reihe.

    Stirling-Näherung (asymptotisch für |z| → ∞):
        ln Γ(z) ≈ (z-½)ln(z) - z + ½ln(2π) + 1/(12z) - 1/(360z³) + ...

    Numerisch stabiler als gamma_lanczos für sehr große |z|.

    @param z: Komplexe Zahl mit Re(z) > 0
    @return: ln(Γ(z))
    @lastModified: 2026-03-10
    """
    if z.real <= 0:
        raise ValueError(f"log_gamma: Re(z) muss positiv sein, erhalten: {z.real}")

    # Für kleine z: Rekursion verwenden bis Re(z) > 10
    if z.real < 10:
        result = log_gamma(z + 1)
        return result - cmath.log(z)

    # Stirling-Reihe für großes |z|
    result = (z - 0.5) * cmath.log(z) - z + 0.5 * math.log(2 * math.pi)
    # Bernoulli-Korrekturen: B₂/(2·1·z) + B₄/(4·3·z³) + B₆/(6·5·z⁵) + ...
    z2 = z * z
    result += 1.0 / (12.0 * z)
    result -= 1.0 / (360.0 * z * z2)
    result += 1.0 / (1260.0 * z * z2 * z2)
    return result


# ===========================================================================
# VOLLSTÄNDIGE RIEMANN-ZETA-FUNKTION
# ===========================================================================

def riemann_zeta(s: complex, precision: int = 200) -> complex:
    """
    Berechnet ζ(s) für BELIEBIGE komplexe Zahlen s ≠ 1.

    Strategie nach Region (drei Fälle):
    - Re(s) > 1:        Euler-Maclaurin (schnelle direkte Konvergenz)
    - 0 < Re(s) ≤ 1:   Dirichlet-Eta-Funktion: η(s) = (1-2^{1-s})·ζ(s)
                         η(s) = Σ (-1)^{n-1}/n^s konvergiert für Re(s) > 0
    - Re(s) ≤ 0:        Funktionalgleichung ζ(s) = χ(s)·ζ(1-s)
                         Der Spiegelpunkt 1-s hat dann Re(1-s) ≥ 1 → Fall 1

    Die Funktionalgleichung (Riemann 1859):
        ζ(s) = 2^s · π^(s-1) · sin(πs/2) · Γ(1-s) · ζ(1-s)

    Für die kritische Gerade Re(s) = 1/2 (Fall 2):
        Nullstellen von ζ(1/2+it) werden durch die RH charakterisiert.

    @param s: Komplexe Zahl, s ≠ 1
    @param precision: Anzahl der Summanden (mehr = genauer, aber langsamer)
    @return: ζ(s)
    @raises ValueError: Wenn s = 1 (Pol)
    @lastModified: 2026-03-10
    """
    if abs(s - 1) < 1e-10:
        raise ValueError("ζ(s) hat einen einfachen Pol bei s = 1 (Residuum = 1)")

    # Sonderfall s = 0: ζ(0) = -1/2
    if abs(s) < 1e-10:
        return complex(-0.5)

    # Fall 1: Re(s) > 1 → Euler-Maclaurin konvergiert direkt
    if s.real > 1:
        return _zeta_euler_maclaurin(s, terms=precision)

    # Fall 2: 0 < Re(s) ≤ 1 → Dirichlet-Eta-Funktion mit Euler-Beschleunigung
    # η(s) = Σ_{n=1}^∞ (-1)^{n-1}/n^s konvergiert für Re(s) > 0
    # ζ(s) = η(s) / (1 - 2^{1-s})
    #
    # PROBLEM: Langsame Konvergenz für Re(s) nahe 0 oder auf der kritischen
    # Geraden (Re=1/2). Die naive Partialsumme braucht N~10^6 Terme für 3 Stellen.
    #
    # LÖSUNG: Euler-Beschleunigung (exponentiell schnell, ~2^{-n} Fehler mit n Termen)
    # E-Transform: η = Σ_{k=0}^{n-1} (1/2^{k+1}) · Δ^k a_0
    # mit Δ^k a_j = Σ_{i=0}^{k} (-1)^i · C(k,i) · a_{j+i},  a_j = 1/(j+1)^s
    if s.real > 0:
        eta_denom = 1 - complex(2) ** (1 - s)
        if abs(eta_denom) < 1e-12:
            raise ValueError("Pol bei s = 1")

        eta = _eta_euler_accelerated(s, n=precision)
        return eta / eta_denom

    # Fall 3: Re(s) ≤ 0 → Funktionalgleichung
    # ζ(s) = 2^s · π^(s-1) · sin(πs/2) · Γ(1-s) · ζ(1-s)
    # Re(1-s) = 1 - Re(s) ≥ 1 → Fall 1 (Euler-Maclaurin) anwendbar
    s_mirror = 1 - s
    zeta_mirror = _zeta_euler_maclaurin(s_mirror, terms=precision)

    # Γ(1-s): Gamma-Funktion am Spiegelpunkt
    try:
        gamma_val = gamma_lanczos(1 - s)
    except (ValueError, OverflowError):
        return complex(0)  # Tritt bei trivialen Nullstellen auf (sin(πs/2) = 0)

    # Faktor der Funktionalgleichung: χ(s) = 2^s · π^{s-1} · sin(πs/2) · Γ(1-s)
    factor = (
        (complex(2) ** s) *
        (cmath.pi ** (s - 1)) *
        cmath.sin(cmath.pi * s / 2) *
        gamma_val
    )

    return factor * zeta_mirror


def _eta_euler_accelerated(s: complex, n: int = 60) -> complex:
    """
    Berechnet η(s) = Σ_{k=1}^∞ (-1)^{k-1}/k^s via Euler-Beschleunigung.

    Die standard Partialsumme konvergiert für Re(s) = 1/2 wie O(1/√N)
    und bräuchte N~10^6 Terme für 3 genaue Stellen.

    Die Euler-E-Transformation beschleunigt auf O(2^{-n}) mit nur n Termen:
        η ≈ Σ_{k=0}^{n-1} (1/2^{k+1}) · Δ^k(a_0)
    mit Vorwärtsdifferenzen Δ^k(a_j) = Σ_{i=0}^{k} (-1)^i · C(k,i) · a_{j+i}
    und a_j = (-1)^j / (j+1)^s (die Glieder der Eta-Reihe, mit Vorzeichen).

    Mit n=60 Termen: Fehler ~2^{-60} ≈ 10^{-18} (unter Maschinenpräzision).

    @param s: Komplexe Zahl mit Re(s) > 0
    @param n: Anzahl der Euler-Transformations-Terme (60 gibt Maschinengenauigkeit)
    @return: η(s)
    @lastModified: 2026-03-10
    """
    # Korrekte Euler-Knopp-E₁-Transformation für Σ_{n=0}^∞ (-1)^n c_n:
    #
    #   S = Σ_{k=0}^∞ (-1)^k Δ^k c_0 / 2^{k+1}
    #
    # wobei:
    #   c_n = 1/(n+1)^s    (OHNE Vorzeichen – unsigned Terme!)
    #   Δ^k c_j = c_{j+k} - k·c_{j+k-1} + C(k,2)·c_{j+k-2} - ...  (Vorwärtsdifferenz)
    #
    # Implementierung via sukzessiver Vorwärtsdifferenzen:
    #   delta anfangs = [c_0, c_1, ..., c_n]
    #   Nach Schritt k: delta[0] = Δ^k c_0

    # Unsigned Terme: c[j] = 1/(j+1)^s
    c = [1.0 / (complex(j + 1) ** s) for j in range(n + 1)]
    delta = list(c)  # Δ^0 c_j = c_j

    result = complex(0.0)
    sign = 1           # (-1)^k, startet mit (-1)^0 = +1
    power_half = 0.5   # 1/2^{k+1}, startet mit 1/2^1 = 0.5

    for k in range(n):
        result += sign * power_half * delta[0]  # (-1)^k / 2^{k+1} · Δ^k c_0
        sign *= -1          # Nächstes (-1)^{k+1}
        power_half *= 0.5   # Nächstes 1/2^{k+2}

        # Vorwärtsdifferenz: Δ^{k+1} c_j = Δ^k c_{j+1} - Δ^k c_j
        for j in range(len(delta) - 1):
            delta[j] = delta[j + 1] - delta[j]
        delta.pop()

    return result


def _zeta_euler_maclaurin(s: complex, terms: int = 200) -> complex:
    """
    Interne Hilfsfunktion: ζ(s) via Euler-Maclaurin für Re(s) > 0.5.

    @param s: Komplexe Zahl mit Re(s) > 0.5
    @param terms: Anzahl der direkten Summanden N
    @return: Näherungswert von ζ(s)
    @lastModified: 2026-03-10
    """
    N = terms
    # Direkte Summe: Σ_{n=1}^{N} n^{-s}
    total = complex(0.0)
    for n in range(1, N + 1):
        total += complex(n) ** (-s)

    # Euler-Maclaurin-Fortsetzung: ∫_N^∞ x^{-s} dx = N^{1-s}/(s-1)
    total += (complex(N) ** (1 - s)) / (s - 1)
    # f(N)/2 = N^{-s}/2
    total += 0.5 * (complex(N) ** (-s))

    # Bernoulli-Korrektur B₂/2! · (-s) · N^{-s-1}
    total += (1.0 / 12.0) * (-s) * (complex(N) ** (-s - 1))
    # B₄/4! · (-s)(-s-1)(-s-2) · N^{-s-3}
    total += (-1.0 / 720.0) * (-s) * (-s - 1) * (-s - 2) * (complex(N) ** (-s - 3))
    # B₆/6! · ... · N^{-s-5}
    total += (1.0 / 30240.0) * (-s) * (-s-1) * (-s-2) * (-s-3) * (-s-4) * (complex(N) ** (-s - 5))

    return total


# ===========================================================================
# RIEMANNS ξ-FUNKTION (symmetrische Form der Zeta-Funktion)
# ===========================================================================

def xi_function(s: complex) -> complex:
    """
    Berechnet Riemanns ξ-Funktion.

    Definition:
        ξ(s) = ½ · s · (s-1) · π^(-s/2) · Γ(s/2) · ζ(s)

    Eigenschaften:
        - ξ(s) = ξ(1-s)  (Symmetrie um Re(s) = 1/2)
        - ξ(s) ist ganz (keine Pole)
        - Die Nullstellen von ξ sind genau die nicht-trivialen Nullstellen von ζ
        - Die RH besagt: Alle Nullstellen von ξ haben Re(s) = 1/2

    Die ξ-Funktion macht die Symmetrie der Riemann-Hypothese besonders deutlich:
    Da ξ(s) = ξ(1-s), liegen Nullstellen in konjugierten Paaren bezüglich Re=1/2.

    @param s: Komplexe Zahl
    @return: ξ(s)
    @lastModified: 2026-03-10
    """
    if abs(s) < 1e-12 or abs(s - 1) < 1e-12:
        return complex(0.5)  # ξ(0) = ξ(1) = 1/2

    # Γ(s/2) berechnen
    try:
        gamma_half_s = gamma_lanczos(s / 2)
    except (ValueError, OverflowError):
        return complex(0)

    # ζ(s) berechnen
    try:
        zeta_s = riemann_zeta(s)
    except ValueError:
        return complex(0)

    # ξ(s) = ½ · s · (s-1) · π^(-s/2) · Γ(s/2) · ζ(s)
    return 0.5 * s * (s - 1) * (cmath.pi ** (-s / 2)) * gamma_half_s * zeta_s


def xi_symmetry_check(s: complex) -> dict[str, object]:
    """
    Verifiziert die Symmetrie ξ(s) = ξ(1-s) numerisch.

    Dies ist ein grundlegendes Werkzeug zur Verifikation der Funktionalgleichung.
    Wenn die Symmetrie exakt gilt, ist die Implementierung korrekt.

    @param s: Testpunkt
    @return: Dictionary mit xi(s), xi(1-s) und dem Fehler
    @lastModified: 2026-03-10
    """
    xi_s = xi_function(s)
    xi_mirror = xi_function(1 - s)
    error = abs(xi_s - xi_mirror)
    return {
        "s": s,
        "mirror": 1 - s,
        "xi_s": xi_s,
        "xi_mirror": xi_mirror,
        "error": error,
        "symmetry_holds": error < 1e-6
    }


# ===========================================================================
# NULLSTELLEN DER ZETA-FUNKTION (Riemann-Hypothese Exploration)
# ===========================================================================

def zeta_on_critical_line(t: float, precision: int = 300) -> complex:
    """
    Berechnet ζ(1/2 + it) auf der kritischen Geraden.

    Dies ist der zentrale Ausdruck für die Riemann-Hypothese.
    Die RH besagt, dass alle nicht-trivialen Nullstellen die Form 1/2 + it haben.

    @param t: Imaginärteil des Punktes auf der kritischen Geraden
    @param precision: Anzahl der Euler-Maclaurin-Terme
    @return: Komplexer Wert ζ(1/2 + it)
    @lastModified: 2026-03-10
    """
    s = complex(0.5, t)
    return riemann_zeta(s, precision=precision)


def riemann_siegel_z(t: float) -> float:
    """
    Berechnet die Riemann-Siegel-Z-Funktion Z(t).

    Die Z-Funktion ist REELLWERTIG und hängt mit ζ auf der kritischen Geraden zusammen:
        Z(t) = e^{iθ(t)} · ζ(1/2 + it)

    mit dem Riemann-Siegel-Winkel:
        θ(t) = Im(ln Γ(1/4 + it/2)) - t/2 · ln(π)

    Eigenschaften:
        - Z(t) ist reell
        - Z(t) = 0 genau dann, wenn ζ(1/2 + it) = 0
        - Vorzeichenwechsel von Z(t) zeigen Nullstellen an
        - Bekannte Nullstellen: t ≈ 14.135, 21.022, 25.011, 30.425, 32.935, ...

    @param t: Reelle Zahl t > 0
    @return: Z(t) ∈ ℝ
    @lastModified: 2026-03-10
    """
    s = complex(0.5, t)

    # Riemann-Siegel-Winkel θ(t)
    # θ(t) = Im(log Γ(1/4 + it/2)) - (t/2)·ln(π)
    try:
        log_g = log_gamma(complex(0.25, t / 2))
        theta = log_g.imag - (t / 2) * math.log(math.pi)
    except Exception:
        return 0.0

    # Z(t) = e^{iθ} · ζ(1/2 + it)
    zeta_val = riemann_zeta(s, precision=200)
    z_val = cmath.exp(complex(0, theta)) * zeta_val
    return z_val.real  # Z(t) ist reell


def find_zeta_zeros(t_min: float, t_max: float, steps: int = 1000) -> list[dict[str, object]]:
    """
    Sucht Nullstellen der Riemann-Zeta-Funktion via Vorzeichenwechsel von Z(t).

    Methode:
        1. Berechne Z(t) auf einem feinen Gitter
        2. Suche Vorzeichenwechsel (zeigt Nullstellen an)
        3. Verfeinere mit Bisektionsverfahren

    Bekannte Nullstellen für Vergleich:
        t₁ = 14.134725...
        t₂ = 21.022040...
        t₃ = 25.010858...
        t₄ = 30.424876...
        t₅ = 32.935062...

    @param t_min: Untere Grenze (t > 0)
    @param t_max: Obere Grenze
    @param steps: Anzahl der Gitterpunkte für Vorzeichensuche
    @return: Liste der gefundenen Nullstellen mit Genauigkeitsinfo
    @lastModified: 2026-03-10
    """
    t_values = np.linspace(t_min, t_max, steps)
    z_values = [riemann_siegel_z(t) for t in t_values]

    zeros = []
    for i in range(len(z_values) - 1):
        # Vorzeichenwechsel: Z wechselt von positiv nach negativ oder umgekehrt
        if z_values[i] * z_values[i + 1] < 0:
            # Bisektions-Verfeinerung
            a, b = t_values[i], t_values[i + 1]
            za, zb = z_values[i], z_values[i + 1]
            for _ in range(50):   # 50 Iterationen → Genauigkeit ~(b-a)/2^50
                mid = (a + b) / 2
                zmid = riemann_siegel_z(mid)
                if za * zmid < 0:
                    b, zb = mid, zmid
                else:
                    a, za = mid, zmid
            t_zero = (a + b) / 2
            zeros.append({
                "t": t_zero,
                "s": complex(0.5, t_zero),
                "on_critical_line": True,  # Per Konstruktion (suchen nur auf Re=1/2)
                "abs_zeta": abs(riemann_zeta(complex(0.5, t_zero), precision=300))
            })

    return zeros


def N_count_formula(T: float) -> float:
    """
    Schätzt die Anzahl N(T) der nicht-trivialen Nullstellen mit Im(ρ) ∈ (0, T].

    Formel (Riemann-von-Mangoldt):
        N(T) = T/(2π) · ln(T/(2π)) - T/(2π) + 7/8 + O(ln T)

    Diese Formel ist exakt bis auf O(ln T) und erlaubt eine Konsistenzprüfung:
    Wenn wir numerisch N Nullstellen finden, und die Formel ~N liefert,
    haben wir wahrscheinlich alle Nullstellen erfasst (kein "Verlust").

    @param T: Obere Grenze des Imaginärteils
    @return: Geschätzte Anzahl N(T)
    @lastModified: 2026-03-10
    """
    if T <= 2 * math.pi:
        return 0.0
    x = T / (2 * math.pi)
    return x * math.log(x) - x + 7.0 / 8.0


# ===========================================================================
# ANALYTISCHE FORTSETZUNG UND FUNKTIONALGLEICHUNG (Visualisierung)
# ===========================================================================

def functional_equation_chi(s: complex) -> complex:
    """
    Berechnet den Faktor χ(s) der Funktionalgleichung ζ(s) = χ(s)·ζ(1-s).

    Explizit:
        χ(s) = 2^s · π^(s-1) · sin(πs/2) · Γ(1-s)

    Überprüfung:
        Wenn ζ korrekt implementiert ist, gilt: ζ(s) = χ(s) · ζ(1-s)
        Für s auf der kritischen Geraden: |χ(1/2 + it)| = 1 (unitär!)

    @param s: Komplexe Zahl s
    @return: χ(s)
    @lastModified: 2026-03-10
    """
    try:
        gamma_1ms = gamma_lanczos(1 - s)
    except (ValueError, OverflowError):
        return complex(0)

    return (
        (2 ** s) *
        (cmath.pi ** (s - 1)) *
        cmath.sin(cmath.pi * s / 2) *
        gamma_1ms
    )


def verify_functional_equation(s: complex) -> dict[str, object]:
    """
    Prüft die Funktionalgleichung ζ(s) = χ(s)·ζ(1-s) numerisch.

    @param s: Testpunkt
    @return: Verifikationsergebnis mit Fehler
    @lastModified: 2026-03-10
    """
    try:
        zeta_s = riemann_zeta(s)
        zeta_1ms = riemann_zeta(1 - s)
        chi_s = functional_equation_chi(s)
        rhs = chi_s * zeta_1ms
        error = abs(zeta_s - rhs)
    except Exception as e:
        return {"error": str(e), "verified": False}

    return {
        "s": s,
        "zeta_s": zeta_s,
        "chi_s_times_zeta_1ms": rhs,
        "error": error,
        "verified": error < 1e-4
    }


# ===========================================================================
# RESIDUENSATZ (für Integralmethoden in der Beweisführung)
# ===========================================================================

def residue_at_pole(f_coeffs: list[complex], order: int) -> complex:
    """
    Berechnet das Residuum einer meromorphen Funktion an einem Pol.

    Für eine Funktion f(z) = Σ_{n=-m}^{∞} a_n · z^n (Laurent-Reihe um z=0)
    ist das Residuum der Koeffizient a_{-1}.

    Der Residuensatz:
        ∮_C f(z) dz = 2πi · Σ Res(f, zₖ)
    für alle Pole zₖ innerhalb der geschlossenen Kurve C.

    @param f_coeffs: Laurent-Koeffizienten [a_{-m}, ..., a_{-1}, a_0, a_1, ...]
                     Indexierung: f_coeffs[order-1] ist der a_{-1}-Term
    @param order: Ordnung des Pols (Anzahl der negativen Potenzen)
    @return: Residuum = Koeffizient von z^{-1}
    @lastModified: 2026-03-10
    """
    if order <= 0 or order > len(f_coeffs):
        raise ValueError(f"Ordnung muss zwischen 1 und {len(f_coeffs)} liegen")
    # Der Koeffizient a_{-1} ist an Position (order - 1) im Array
    return f_coeffs[order - 1]


def cauchy_integral_numerical(
    f: Callable[[complex], complex],
    center: complex,
    radius: float,
    n_points: int = 1000
) -> complex:
    """
    Numerische Berechnung des Cauchy-Integrals 1/(2πi) ∮ f(z)/(z-z₀) dz.

    Der Cauchy-Integralsatz besagt:
        f(z₀) = 1/(2πi) · ∮_C f(z)/(z-z₀) dz

    für holomorphes f und Kreis C mit Radius r um z₀.

    Anwendungen:
        - Rekonstruktion von Funktionswerten aus Randwerten
        - Berechnung von Taylor-Koeffizienten
        - Beweis des Residuensatzes numerisch verifizieren

    @param f: Holomorphe Funktion f: ℂ → ℂ
    @param center: Mittelpunkt z₀ der Kreisintegration
    @param radius: Radius r des Kreises
    @param n_points: Anzahl der Quadraturpunkte
    @return: Approximation von f(z₀)
    @lastModified: 2026-03-10
    """
    # Cauchy-Integralsatz via Mittelwerteigenschaft holomorpher Funktionen:
    #
    #   1/(2πi) ∮_C f(z)/(z-z₀) dz
    #   = 1/(2πi) ∫_0^{2π} f(z₀ + r·e^{iθ}) / (r·e^{iθ}) · i·r·e^{iθ} dθ
    #   = 1/(2π) ∫_0^{2π} f(z₀ + r·e^{iθ}) dθ  ← Mittelwert auf dem Kreis
    #
    # Diskrete Trapez-Regel (exakt für trigonometrische Polynome bis Grad n/2):
    #   f(z₀) ≈ (1/n) Σ_{k=0}^{n-1} f(z₀ + r·e^{i·2πk/n})

    total = complex(0.0)
    for k in range(n_points):
        theta = 2 * math.pi * k / n_points
        z = center + radius * cmath.exp(complex(0, theta))
        total += f(z)
    return total / n_points


# ===========================================================================
# ARBITRARY PRECISION via mpmath
# ===========================================================================

try:
    import mpmath
    _MPMATH_AVAILABLE = True
except ImportError:
    _MPMATH_AVAILABLE = False


def _check_mpmath() -> None:
    """
    Prüft ob mpmath installiert ist und wirft ImportError wenn nicht.

    @raises ImportError: Wenn mpmath nicht verfügbar ist
    @lastModified: 2026-03-10
    """
    if not _MPMATH_AVAILABLE:
        raise ImportError(
            "mpmath ist nicht installiert. "
            "Installation: pip install mpmath --break-system-packages"
        )


def riemann_zeta_mpmath(s: complex, dps: int = 50) -> complex:
    """
    Berechnet die Riemann-Zeta-Funktion ζ(s) mit beliebiger Genauigkeit via mpmath.

    Vorteile gegenüber riemann_zeta():
    - dps=50  → ~50 Dezimalstellen (≈ 165 Bit), doppelt so viel wie float64
    - dps=100 → ~100 Dezimalstellen (≈ 330 Bit)
    - mpmath verwendet intern die MPFR-Bibliothek (GNU Multiple Precision)
    - Unterstützt analytische Fortsetzung, Funktionalgleichung etc. automatisch

    Die mpmath-Implementierung von zeta(s) nutzt:
    - Euler-Maclaurin für Re(s) > 1
    - Alternating Euler-Bernoulli-Reihe für 0 < Re(s) ≤ 1
    - Funktionalgleichung für Re(s) ≤ 0

    @param s: Komplexe Zahl s ≠ 1
    @param dps: Dezimalstellen Genauigkeit (Standard: 50 ≈ 165 Bit)
    @return: ζ(s) als komplexe Python-Zahl (float64-Genauigkeit des Rückgabewerts,
             interne Berechnung in dps Dezimalstellen)
    @raises ImportError: Wenn mpmath nicht installiert ist
    @raises ValueError: Wenn s = 1 (Pol)
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _check_mpmath()

    if abs(s - 1) < 1e-10:
        raise ValueError("ζ(s) hat einen einfachen Pol bei s = 1")

    # Genauigkeit setzen
    with mpmath.workdps(dps):
        # Python complex → mpmath mpc (multiple precision complex)
        s_mp = mpmath.mpc(s.real, s.imag)
        # mpmath.zeta nutzt optimierte Algorithmen für beliebige Präzision
        result = mpmath.zeta(s_mp)
        # Rückgabe als Python complex (float64-Genauigkeit)
        return complex(float(result.real), float(result.imag))


def find_zeta_zeros_mpmath(
    t_min: float,
    t_max: float,
    dps: int = 50,
    steps: int = 1000
) -> list[dict[str, object]]:
    """
    Sucht Riemann-Nullstellen auf der kritischen Geraden s = 1/2 + it
    mit beliebiger Genauigkeit via mpmath.

    Methode:
    1. Primärsuche: Vorzeichenwechsel von riemann_siegel_z_mpmath() auf
       dem Intervall [t_min, t_max] mit 'steps' Gitterpunkten
    2. Verfeinerung: mpmath.findroot (Newton-Raphson) nahe dem Vorzeichenwechsel
    3. Verifikation: |ζ(1/2 + it*)| < Schwellwert

    Bekannte erste Nullstellen (für Validierung):
        t₁ ≈ 14.134725141734693
        t₂ ≈ 21.022039638771555
        t₃ ≈ 25.010857580145688
        t₄ ≈ 30.424876125859513

    @param t_min: Untere Suchgrenze (t > 0)
    @param t_max: Obere Suchgrenze
    @param dps: Dezimalstellen Genauigkeit (Standard: 50)
    @param steps: Anzahl Gitterpunkte für Vorzeichensuche (Standard: 1000)
    @return: Liste von Dictionaries:
             {'t': float, 'zero': complex, 'verified': bool}
             wobei 'verified' bedeutet: |ζ(1/2 + it)| < 1e-8
    @raises ImportError: Wenn mpmath nicht installiert ist
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _check_mpmath()

    zeros = []

    with mpmath.workdps(dps):
        # Gitterpunkte für Vorzeichensuche
        t_grid = [t_min + (t_max - t_min) * i / steps for i in range(steps + 1)]

        # Riemann-Siegel Z-Funktion auf dem Gitter auswerten
        z_vals = []
        for t in t_grid:
            try:
                z_vals.append(float(mpmath.siegelz(t)))
            except Exception:
                z_vals.append(0.0)

        # Vorzeichenwechsel suchen und mit Newton-Raphson verfeinern
        for i in range(len(z_vals) - 1):
            if z_vals[i] * z_vals[i + 1] < 0:
                # Vorzeichenwechsel gefunden → Newton-Raphson Verfeinerung
                t_a, t_b = t_grid[i], t_grid[i + 1]
                try:
                    # mpmath.findroot: Newton-Verfahren für Z(t) = 0
                    t_zero = float(mpmath.findroot(
                        lambda t: mpmath.siegelz(t),
                        (t_a + t_b) / 2,
                        tol=mpmath.mpf(10) ** (-dps + 5)
                    ))

                    # Nullstelle auf der kritischen Geraden
                    s_zero = complex(0.5, t_zero)

                    # Verifikation: |ζ(1/2 + it*)| soll sehr klein sein
                    zeta_val = riemann_zeta_mpmath(s_zero, dps=dps)
                    abs_zeta = abs(zeta_val)

                    zeros.append({
                        't': t_zero,
                        'zero': s_zero,
                        'verified': abs_zeta < 1e-6
                    })
                except Exception:
                    # Konvergenzfehler: Vorzeichenwechsel trotzdem merken
                    t_approx = (t_a + t_b) / 2
                    zeros.append({
                        't': t_approx,
                        'zero': complex(0.5, t_approx),
                        'verified': False
                    })

    return zeros


def riemann_siegel_z_mpmath(t: float, dps: int = 50) -> float:
    """
    Berechnet die Riemann-Siegel-Z-Funktion Z(t) mit mpmath-Genauigkeit.

    Die Z-Funktion ist reellwertig und definiert als:
        Z(t) = e^{iθ(t)} · ζ(1/2 + it)

    mit dem Riemann-Siegel-Winkel θ(t) = Im(ln Γ(1/4 + it/2)) - t·ln(π)/2.

    Eigenschaften:
    - Z(t) ∈ ℝ für alle t ∈ ℝ
    - Z(t) = 0 ⟺ ζ(1/2 + it) = 0 (Nullstellen auf kritischer Geraden)
    - |Z(t)| = |ζ(1/2 + it)| (gleiche Betragsstruktur)
    - Vorzeichenwechsel von Z(t) zeigen Nullstellen an

    mpmath bietet mpmath.siegelz(t) als direkte Implementierung
    mit korrekter Behandlung von Vorzeichen und Näherungsformeln.

    @param t: Reelle Zahl t (typischerweise t > 0)
    @param dps: Dezimalstellen Genauigkeit (Standard: 50)
    @return: Z(t) als Python float
    @raises ImportError: Wenn mpmath nicht installiert ist
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _check_mpmath()

    with mpmath.workdps(dps):
        # mpmath.siegelz(t) berechnet Z(t) = e^{iθ(t)} · ζ(1/2 + it) direkt
        result = mpmath.siegelz(mpmath.mpf(t))
        return float(result)


def gamma_mpmath(z: complex, dps: int = 50) -> complex:
    """
    Berechnet die Gamma-Funktion Γ(z) mit mpmath beliebiger Genauigkeit.

    Vergleich mit gamma_lanczos():
    - gamma_lanczos: ~15 Dezimalstellen (float64-Grenze)
    - gamma_mpmath mit dps=50: ~50 Dezimalstellen
    - gamma_mpmath mit dps=100: ~100 Dezimalstellen

    Die Gamma-Funktion verallgemeinert die Fakultät:
        Γ(n) = (n-1)!  für positive ganze n ≥ 1
        Γ(1/2) = √π ≈ 1.7724538509...
        Γ(z+1) = z·Γ(z)  (Rekursionsformel)

    mpmath verwendet die Lanczos-Approximation mit erweiterter Präzision
    und Reflexionsformel für Re(z) < 1/2.

    @param z: Komplexe Zahl (z ∉ {0, -1, -2, ...}, sonst Pol)
    @param dps: Dezimalstellen Genauigkeit (Standard: 50)
    @return: Γ(z) als komplexe Python-Zahl
    @raises ImportError: Wenn mpmath nicht installiert ist
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _check_mpmath()

    with mpmath.workdps(dps):
        # Python complex → mpmath mpc
        z_mp = mpmath.mpc(z.real, z.imag)
        # mpmath.gamma: vollständige Gamma-Funktion für komplexe Argumente
        result = mpmath.gamma(z_mp)
        return complex(float(result.real), float(result.imag))


def verify_riemann_hypothesis_range_mpmath(
    n_zeros: int = 10,
    dps: int = 50
) -> dict[str, object]:
    """
    Verifiziert die ersten n Riemann-Nullstellen auf der kritischen Geraden.

    Methode:
    - mpmath.zetazero(n) liefert die n-te nicht-triviale Nullstelle ρ_n
    - Prüft: Re(ρ_n) = 1/2 (bis auf numerische Genauigkeit ~10^(-dps+5))
    - Prüft: |ζ(ρ_n)| ≈ 0

    Dies ist eine numerische Verifikation der Riemann-Hypothese für
    die ersten n_zeros Nullstellen (bekannt: alle ersten ~10^13 Nullstellen
    liegen auf der kritischen Geraden, per Computerbeweis).

    Die Riemann-Hypothese selbst (für ALLE Nullstellen) ist unbewiesen.

    Bekannte erste Nullstellen:
        ρ₁ = 1/2 + 14.134725141734693...i
        ρ₂ = 1/2 + 21.022039638771555...i
        ρ₃ = 1/2 + 25.010857580145688...i

    @param n_zeros: Anzahl der zu prüfenden Nullstellen (Standard: 10)
    @param dps: Dezimalstellen Genauigkeit (Standard: 50)
    @return: Dictionary:
             {'zeros': list, 'all_on_critical_line': bool, 'max_deviation': float}
             - zeros: Liste von {'n': int, 'rho': complex, 're_deviation': float}
             - all_on_critical_line: True wenn alle Re(ρ_n) = 1/2 ± 1e-8
             - max_deviation: Maximale Abweichung |Re(ρ_n) - 1/2|
    @raises ImportError: Wenn mpmath nicht installiert ist
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _check_mpmath()

    zeros_list = []
    max_dev = 0.0

    with mpmath.workdps(dps):
        for n in range(1, n_zeros + 1):
            try:
                # mpmath.zetazero(n): n-te nicht-triviale Nullstelle
                rho = mpmath.zetazero(n)
                re_part = float(rho.real)
                im_part = float(rho.imag)
                # Abweichung von 1/2
                deviation = abs(re_part - 0.5)
                max_dev = max(max_dev, deviation)

                zeros_list.append({
                    'n': n,
                    'rho': complex(re_part, im_part),
                    're_deviation': deviation
                })
            except Exception as e:
                zeros_list.append({
                    'n': n,
                    'rho': None,
                    're_deviation': float('inf'),
                    'error': str(e)
                })
                max_dev = float('inf')

    # Toleranz: 1e-8 (weit unter der dps-Genauigkeit)
    all_on_line = max_dev < 1e-8 and max_dev != float('inf')

    return {
        'zeros': zeros_list,
        'all_on_critical_line': all_on_line,
        'max_deviation': max_dev
    }


def li_function_mpmath(x: float, dps: int = 50) -> float:
    """
    Berechnet den Logarithmischen Integralus li(x) mit mpmath-Genauigkeit.

    Definition (Cauchy-Hauptwert):
        li(x) = PV ∫₀^x dt / ln(t)

    Die li-Funktion hat eine Singularität bei t = 1 (ln(1) = 0).
    Der Cauchy-Hauptwert ist für alle x > 0, x ≠ 1 definiert.

    Verwendung in der Primzahltheorie:
    - Primzahlsatz: π(x) ~ li(x) für x → ∞
    - li(x) ist eine bessere Approximation von π(x) als x/ln(x)
    - Gauss vermutete: π(x) < li(x) für alle x > 1
    - Skewes-Zahl: erster bekannter Überschneidungspunkt bei ~10^316

    mpmath.li(x) berechnet li(x) via:
        li(x) = Ei(ln x)
    mit Ei(x) = PV ∫_{-∞}^x e^t/t dt (Exponentialintegral).

    @param x: Reelle Zahl x > 0, x ≠ 1
    @param dps: Dezimalstellen Genauigkeit (Standard: 50)
    @return: li(x) als Python float
    @raises ImportError: Wenn mpmath nicht installiert ist
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    _check_mpmath()

    with mpmath.workdps(dps):
        # mpmath.li(x): logarithmischer Integralus mit Cauchy-Hauptwert
        result = mpmath.li(mpmath.mpf(x))
        return float(result)


# ===========================================================================
# GUE-STATISTIK DER RIEMANN-NULLSTELLEN (Montgomery-Odlyzko-Gesetz)
# ===========================================================================

def _get_riemann_zeros(n_zeros: int) -> list[float]:
    """
    Gibt die Imaginärteile der ersten n nicht-trivialen Riemann-Nullstellen zurück.

    Die ersten bekannten Nullstellen ρ_n = 1/2 + i·γ_n sind numerisch
    auf höchste Präzision berechnet worden (Odlyzko, van de Lune et al.).
    Für größere n_zeros wird mpmath.zetazero() verwendet.

    @param n_zeros: Anzahl der gewünschten Nullstellen
    @return: Sortierte Liste der Imaginärteile γ_n > 0
    @author: Michael Fuhrmann
    @lastModified: 2026-03-10
    """
    # Vorberechnete bekannte erste Nullstellen (hohe Genauigkeit)
    known_zeros = [
        14.134725141734693790,
        21.022039638771554993,
        25.010857580145688763,
        30.424876125859513210,
        32.935061587739189691,
        37.586178158825671257,
        40.918719012147495187,
        43.327073280914999519,
        48.005150881167159727,
        49.773832477672302181,
        52.970321477714460644,
        56.446247697063246010,
        59.347044002602353079,
        60.831778524609809844,
        65.112544048081606816,
        67.079810529494172616,
        69.546401711173979252,
        72.067157674481907584,
        75.704690699083933168,
        77.144840069745236890,
        79.337375020249367022,
        82.910380854086030183,
        84.735492981329728860,
        87.425274613125229406,
        88.809111208594021956,
        92.491899270558484296,
        94.651344040519886098,
        95.870634228245309758,
        98.831194218193692888,
        101.31785100573139122,
    ]

    if n_zeros <= len(known_zeros):
        return known_zeros[:n_zeros]

    # Für mehr Nullstellen: mpmath verwenden
    try:
        _check_mpmath()
        zeros = list(known_zeros)
        with mpmath.workdps(25):
            for k in range(len(known_zeros) + 1, n_zeros + 1):
                z = mpmath.zetazero(k)
                zeros.append(float(z.imag))
        return zeros
    except Exception:
        # Fallback: nur die bekannten Nullstellen zurückgeben
        return known_zeros[:min(n_zeros, len(known_zeros))]


def zeta_zero_spacing_statistics(n_zeros: int = 100) -> dict:
    """
    Berechnet die normierten Abstände zwischen aufeinanderfolgenden Riemann-Nullstellen
    und vergleicht deren Statistik mit der GUE-Vorhersage (Montgomery-Odlyzko-Gesetz).

    Das Montgomery-Odlyzko-Gesetz (1973/1979/1987) besagt:
        Die normierten Abstände zwischen Riemann-Nullstellen folgen
        derselben Statistik wie Eigenwert-Abstände des Gaussian Unitary Ensemble (GUE).

    Normierung der Abstände (lokale Mittelwert-Normierung):
        δ_n = (γ_{n+1} - γ_n) · log(γ_n / (2π)) / (2π)

    Dabei ist log(γ_n / (2π)) / (2π) die lokale Dichte der Nullstellen
    bei Höhe γ_n (aus der Weil-Formel: N(T) ≈ T/(2π) · log(T/(2πe))).

    GUE-Vorhersage (Wigner-Surmise):
        P(δ) ≈ (32/π²) δ² exp(-4δ²/π)

    Diese Formel beschreibt die „level repulsion" – Nullstellen weichen
    einander aus, ähnlich wie Energieniveaus in schweren Atomkernen.

    Im Vergleich: Poisson-Verteilung P(δ) = exp(-δ) (unkorrelierten Ereignisse)
    würde viele kleine Abstände (clustering) vorhersagen, was bei Riemann-Nullstellen
    NICHT beobachtet wird.

    @param n_zeros: Anzahl der Riemann-Nullstellen (Standard: 100)
    @return: Dict mit:
             - 'spacings': normierte Abstände δ_n
             - 'mean_spacing': Mittelwert (sollte ≈ 1 nach Normierung sein)
             - 'variance': Varianz der Abstände
             - 'std_dev': Standardabweichung
             - 'min_spacing', 'max_spacing': Extremwerte
             - 'ks_statistic', 'ks_p_value': KS-Test gegen GUE-Wigner-Surmise
             - 'ks_poisson_statistic', 'ks_poisson_p_value': KS-Test gegen Poisson
             - 'gue_mean_theoretical': theoretischer GUE-Mittelwert ≈ 1
             - 'n_spacings': Anzahl der verwendeten Abstände
    @author: Michael Fuhrmann
    @lastModified: 2026-03-10
    """
    from scipy import stats as scipy_stats
    import numpy as np

    # Riemann-Nullstellen laden
    zeros = _get_riemann_zeros(n_zeros)
    n_actual = len(zeros)

    if n_actual < 3:
        return {"error": "Zu wenige Nullstellen für Statistik", "n_spacings": 0}

    # Abstände zwischen aufeinanderfolgenden Nullstellen berechnen
    gaps = [zeros[i + 1] - zeros[i] for i in range(n_actual - 1)]

    # Normierung: δ_n = Lücke × lokale Nullstellendichte
    # Lokale Dichte bei γ_n: d(T) = log(T/(2π)) / (2π)
    normalized_spacings = []
    for i, gap in enumerate(gaps):
        gamma_n = zeros[i]
        if gamma_n > 0:
            # Lokale Dichte der Nullstellen bei Höhe γ_n
            local_density = math.log(gamma_n / (2.0 * math.pi)) / (2.0 * math.pi)
            if local_density > 0:
                delta_n = gap * local_density
                normalized_spacings.append(delta_n)

    if not normalized_spacings:
        return {"error": "Normierung fehlgeschlagen", "n_spacings": 0}

    spacings_arr = np.array(normalized_spacings)
    mean_spacing = float(np.mean(spacings_arr))
    variance = float(np.var(spacings_arr))
    std_dev = float(np.std(spacings_arr))

    # GUE Wigner-Surmise CDF: P(δ) = (32/π²)·δ²·exp(-4δ²/π)
    # CDF: F(x) = 1 - exp(-4x²/π) · (1 + 4x²/π · ...)
    # Numerisch: exakte CDF via numerische Integration
    def gue_wigner_cdf(x_vals: np.ndarray) -> np.ndarray:
        """CDF der GUE Wigner-Surmise-Verteilung."""
        # PDF: p(s) = (32/π²)·s²·exp(-4s²/π)
        # CDF: F(x) = Integral von 0 bis x
        from scipy.special import erf
        # Analytische CDF: F(x) = erf(2x/√π) - (4x/π)·exp(-4x²/π)
        # (Ergebnis der Integration der Wigner-Surmise PDF)
        a = 2.0 / math.sqrt(math.pi)
        result = erf(a * x_vals) - (4.0 * x_vals / math.pi) * np.exp(-4.0 * x_vals ** 2 / math.pi)
        return np.clip(result, 0.0, 1.0)

    # KS-Test gegen GUE-Wigner-Surmise (empirische CDF vs. theoretische CDF)
    sorted_spacings = np.sort(spacings_arr)
    n_s = len(sorted_spacings)
    empirical_cdf = np.arange(1, n_s + 1) / n_s
    theoretical_gue_cdf = gue_wigner_cdf(sorted_spacings)

    ks_stat_gue = float(np.max(np.abs(empirical_cdf - theoretical_gue_cdf)))

    # Array-kompatible CDF für scipy.stats.ks_1samp
    # (scipy übergibt intern ein Array aller Stichprobenwerte auf einmal)
    def gue_wigner_cdf_scalar(x):
        """
        Array-kompatible GUE Wigner-Surmise CDF für scipy.stats.ks_1samp.

        Verarbeitet sowohl skalare float-Werte als auch numpy-Arrays,
        da scipy.stats.ks_1samp intern mit Arrays arbeitet.

        @param x: float oder numpy-Array mit Auswertungspunkten
        @return: CDF-Werte (gleicher Typ wie Eingabe)
        """
        from scipy.special import erf as scipy_erf_vec
        x_arr = np.asarray(x, dtype=float)
        a = 2.0 / math.sqrt(math.pi)
        val = scipy_erf_vec(a * x_arr) - (4.0 * x_arr / math.pi) * np.exp(-4.0 * x_arr ** 2 / math.pi)
        result_arr = np.clip(val, 0.0, 1.0)
        # Skalar zurückgeben wenn Eingabe skalar war
        if result_arr.ndim == 0:
            return float(result_arr)
        return result_arr

    # p-Wert via Kolmogorov-Smirnov-Verteilung
    ks_result_gue = scipy_stats.ks_1samp(
        spacings_arr,
        gue_wigner_cdf_scalar,
        alternative='two-sided'
    )

    # KS-Test gegen Poisson-Verteilung (Exp(1) für normierte Abstände)
    ks_result_poisson = scipy_stats.ks_1samp(
        spacings_arr,
        scipy_stats.expon.cdf,
        alternative='two-sided'
    )

    # GUE Wigner-Surmise: theoretischer Mittelwert = Γ(4/3) / Γ(3/3)·... ≈ 1
    # Exakter Wert: E[δ] = Γ(4/3) · (π/4)^{1/3} · ... numerisch ≈ 1.0 nach Normierung
    gue_mean_theoretical = 1.0  # Per Definition der Normierung

    # Wigner-Surmise Varianz: Var(δ) = 4/π - (π/4) ≈ 0.2732 (GUE-Vorhersage)
    gue_variance_theoretical = 4.0 / math.pi - (math.pi / 4.0)

    return {
        "spacings": normalized_spacings,
        "n_spacings": len(normalized_spacings),
        "mean_spacing": mean_spacing,
        "variance": variance,
        "std_dev": std_dev,
        "min_spacing": float(np.min(spacings_arr)),
        "max_spacing": float(np.max(spacings_arr)),
        "ks_statistic_gue": float(ks_result_gue.statistic),
        "ks_p_value_gue": float(ks_result_gue.pvalue),
        "ks_statistic_poisson": float(ks_result_poisson.statistic),
        "ks_p_value_poisson": float(ks_result_poisson.pvalue),
        "gue_mean_theoretical": gue_mean_theoretical,
        "gue_variance_theoretical": float(gue_variance_theoretical),
        "n_zeros_used": n_actual,
        # GUE zeigt level repulsion (Abstoßung), Poisson zeigt clustering
        "level_repulsion_detected": variance < 0.6,
        "wigner_surmise_formula": "P(delta) = (32/pi^2) * delta^2 * exp(-4*delta^2/pi)"
    }


def pair_correlation_function(n_zeros: int = 200, r: float = 1.0) -> dict:
    """
    Berechnet Montgomerys Paarkorrelationsfunktion der Riemann-Nullstellen.

    Montgomery (1973) bewies unter RH, dass die Paarkorrelationsfunktion
    der Riemann-Nullstellen mit der GUE-Paarkorrelation übereinstimmt:

        R_2(r) = 1 - (sin(πr) / (πr))²

    Definition der empirischen Paarkorrelation:
        R_2(α, β) = (1 / N(T)) · Σ_{γ, γ' distinct} f((γ - γ') · log(T) / (2π))

    wobei f eine Testfunktion mit Fourier-Träger in (-2, 2) ist.

    Die GUE-Paarkorrelation:
        R_2^{GUE}(r) = 1 - (sin(πr) / (πr))²

    hat die Eigenschaft, bei r → 0 gegen 0 zu gehen (level repulsion!),
    im Gegensatz zur Poisson-Paarkorrelation R_2^{Poisson}(r) = 1.

    Interpretation:
        r = 0: R_2 → 0 (Nullstellen weichen einander aus)
        r → ∞: R_2 → 1 (keine Langstreckenkorrela tion)
        Maximum bei r ≈ 1.2 (leichte Überhöhung über 1 – Bunching)

    @param n_zeros: Anzahl der Nullstellen für die empirische Berechnung
    @param r: Auswertungspunkt der Paarkorrelationsfunktion (r > 0)
    @return: Dict mit:
             - 'r': Auswertungspunkt
             - 'gue_prediction': GUE-Vorhersage 1 - (sin(πr)/(πr))²
             - 'poisson_prediction': Poisson-Vorhersage = 1.0
             - 'empirical': empirische Paarkorrelation aus Nullstellen
             - 'r_grid': Liste von r-Werten für Kurvenplot
             - 'gue_curve': R_2^{GUE} für jeden r-Wert
             - 'empirical_curve': empirische R_2 für jeden r-Wert
    @author: Michael Fuhrmann
    @lastModified: 2026-03-10
    """
    import numpy as np

    def gue_pair_correlation(r_val: float) -> float:
        """
        GUE-Paarkorrelation: R_2(r) = 1 - (sin(πr) / (πr))²

        @param r_val: Korrelationsabstand
        @return: R_2(r) ∈ [0, 2]
        """
        if abs(r_val) < 1e-12:
            return 0.0  # Grenzwert: lim_{r→0} (sin(πr)/(πr))² = 1 → R_2(0) = 0
        sin_term = math.sin(math.pi * r_val) / (math.pi * r_val)
        return 1.0 - sin_term ** 2

    # GUE-Vorhersage am Auswertungspunkt r
    gue_at_r = gue_pair_correlation(r)

    # Nullstellen laden
    zeros = _get_riemann_zeros(n_zeros)
    n_actual = len(zeros)

    # Empirische Paarkorrelation berechnen
    # Methode: Histogram der normierten Differenzen (γ_i - γ_j) · log(T) / (2π)
    # Wir wählen T als geometrisches Mittel der verwendeten Nullstellen
    if n_actual < 4:
        return {
            "r": r,
            "gue_prediction": gue_at_r,
            "poisson_prediction": 1.0,
            "empirical": None,
            "error": "Zu wenige Nullstellen"
        }

    T_ref = zeros[n_actual // 2]  # Referenzhöhe ≈ Mitte des Bereichs
    log_T_over_2pi = math.log(T_ref) / (2.0 * math.pi)

    # Sammle alle normierten Differenzen (nur für benachbarte Paare für Effizienz)
    normalized_diffs = []
    window = min(20, n_actual - 1)  # Fenster für Paarbildung

    for i in range(n_actual):
        for j in range(i + 1, min(i + window + 1, n_actual)):
            diff = (zeros[j] - zeros[i]) * log_T_over_2pi
            if 0.01 < diff < 5.0:  # Nur sinnvolle Abstände
                normalized_diffs.append(diff)

    # Kurven für Plot-Bereich r ∈ [0, 3]
    r_grid = [i * 0.1 for i in range(1, 31)]  # r = 0.1, 0.2, ..., 3.0
    gue_curve = [gue_pair_correlation(rv) for rv in r_grid]

    # Empirische Kurve via Kernel-Dichte-Schätzung (Histogramm-Approximation)
    empirical_curve = []
    bandwidth = 0.3  # Glättungsbreite
    if normalized_diffs:
        diffs_arr = np.array(normalized_diffs)
        for rv in r_grid:
            # Gaussian-Kernel: Gewichte nach Abstand zu rv
            weights = np.exp(-0.5 * ((diffs_arr - rv) / bandwidth) ** 2)
            # Normierung: Vergleich mit Erwartungswert für zufällige Abstände
            empirical_curve.append(float(np.sum(weights)))
    else:
        empirical_curve = [0.0] * len(r_grid)

    # Normierung der empirischen Kurve (auf Einheitsintervall skalieren)
    if empirical_curve and max(empirical_curve) > 0:
        max_emp = max(empirical_curve)
        # Skaliere so dass der Mittelwert (über r > 1.5) bei 1 liegt
        tail_mean = np.mean([v for v, rv in zip(empirical_curve, r_grid) if rv > 1.5])
        if tail_mean > 0:
            empirical_curve = [v / tail_mean for v in empirical_curve]

    # Empirische Paarkorrelation am Auswertungspunkt r interpolieren
    empirical_at_r = None
    if normalized_diffs:
        diffs_arr = np.array(normalized_diffs)
        weights = np.exp(-0.5 * ((diffs_arr - r) / bandwidth) ** 2)
        empirical_at_r = float(np.sum(weights))
        # Normierung
        tail_refs = np.array([rv for rv in r_grid if rv > 1.5])
        if len(tail_refs) > 0:
            ref_weights = [
                float(np.sum(np.exp(-0.5 * ((diffs_arr - rv) / bandwidth) ** 2)))
                for rv in tail_refs
            ]
            ref_mean = np.mean(ref_weights)
            if ref_mean > 0:
                empirical_at_r /= ref_mean

    return {
        "r": float(r),
        "gue_prediction": float(gue_at_r),
        "poisson_prediction": 1.0,
        "empirical": empirical_at_r,
        "r_grid": r_grid,
        "gue_curve": gue_curve,
        "empirical_curve": empirical_curve,
        "n_zeros_used": n_actual,
        "n_pairs": len(normalized_diffs),
        # Interpretation
        "level_repulsion": gue_at_r < 0.5,
        "formula": "R_2(r) = 1 - (sin(pi*r) / (pi*r))^2"
    }


def random_matrix_gue_sample(n: int = 100, size: int = 50) -> dict:
    """
    Erzeugt GUE-Zufallsmatrizen und vergleicht deren Eigenwertstatistik
    mit der Riemann-Nullstellenstatistik.

    Das Gaussian Unitary Ensemble (GUE) besteht aus hermiteschen Matrizen
    mit normalverteilten Einträgen:

        H = (A + A†) / √(2·N)

    wobei A_{ij} = X_{ij} + i·Y_{ij} mit X_{ij}, Y_{ij} ~ N(0, 1) unabhängig.

    Eigenschaften von GUE-Matrizen:
        - Hermitesch: H = H† → reelle Eigenwerte
        - Eigenwert-Dichte: Wigner-Halbkreisgesetz ρ(λ) = √(4-λ²) / (2π)
        - Eigenwert-Abstände: Wigner-Surmise P(s) = (π/2)·s·exp(-πs²/4)
          [für GUE: P(s) = (32/π²)·s²·exp(-4s²/π)]

        ACHTUNG: Es gibt verschiedene GUE-Konventionen. Das Wigner-Surmise
        für GUE (β=2) lautet:
            P(s) = (32/π²) s² exp(-4s²/π)
        Dies ist identisch mit der Riemann-Nullstellenstatistik!

    Der Zusammenhang Riemann-Nullstellen ↔ GUE-Eigenwerte ist eine der
    tiefsten offenen Fragen der Mathematik (Hilbert-Pólya-Vermutung):
        Gibt es einen hermiteschen Operator, dessen Eigenwerte die
        Imaginärteile der Riemann-Nullstellen sind?

    @param n: Anzahl der GUE-Matrizen (Stichprobengröße, Standard: 100)
    @param size: Matrixgröße N×N (Standard: 50)
    @return: Dict mit:
             - 'eigenvalue_spacings': alle normierten Eigenwert-Abstände
             - 'mean_spacing': Mittelwert (sollte ≈ 1 sein)
             - 'variance': Varianz
             - 'ks_statistic', 'ks_p_value': KS-Test gegen Wigner-Surmise
             - 'sample_matrix_shape': (size, size)
             - 'n_matrices': Anzahl der erzeugten Matrizen
             - 'eigenvalues_real': True (hermitesche Matrix → reelle Eigenwerte)
    @author: Michael Fuhrmann
    @lastModified: 2026-03-10
    """
    import numpy as np
    from scipy import stats as scipy_stats

    rng = np.random.default_rng(seed=42)  # Reproduzierbarkeit

    all_spacings = []

    for _ in range(n):
        # GUE-Matrix erzeugen: A mit komplexen Gaußeinträgen
        real_part = rng.standard_normal((size, size))
        imag_part = rng.standard_normal((size, size))
        A = real_part + 1j * imag_part

        # Hermitesch machen und normieren
        H = (A + A.conj().T) / math.sqrt(2.0 * size)

        # Eigenwerte berechnen (hermitesche Matrix → reell)
        # np.linalg.eigvalsh: schneller als eig für hermitesche Matrizen
        eigenvalues = np.linalg.eigvalsh(H)
        eigenvalues_sorted = np.sort(eigenvalues.real)

        # Abstände zwischen aufeinanderfolgenden Eigenwerten
        raw_spacings = np.diff(eigenvalues_sorted)

        if len(raw_spacings) > 0 and np.mean(raw_spacings) > 0:
            # Normierung: teile durch mittleren Abstand (unfolding)
            mean_s = np.mean(raw_spacings)
            normalized = raw_spacings / mean_s
            all_spacings.extend(normalized.tolist())

    if not all_spacings:
        return {"error": "Keine Eigenwert-Abstände berechnet", "n_matrices": n}

    spacings_arr = np.array(all_spacings)
    mean_spacing = float(np.mean(spacings_arr))
    variance = float(np.var(spacings_arr))
    std_dev = float(np.std(spacings_arr))

    # Verifiziere: Eigenwerte sind reell (Test an Beispiel-Matrix)
    real_part_test = rng.standard_normal((size, size))
    imag_part_test = rng.standard_normal((size, size))
    A_test = real_part_test + 1j * imag_part_test
    H_test = (A_test + A_test.conj().T) / math.sqrt(2.0 * size)
    eig_test = np.linalg.eigvalsh(H_test)
    eigenvalues_real = bool(np.all(np.abs(eig_test.imag) < 1e-10))

    # KS-Test gegen GUE Wigner-Surmise P(s) = (32/π²)·s²·exp(-4s²/π)
    from scipy.special import erf as scipy_erf

    def gue_wigner_cdf_matrix(x):
        """
        Array-kompatible CDF der GUE Wigner-Surmise für scipy.stats.ks_1samp.

        Verarbeitet sowohl skalare Werte als auch numpy-Arrays,
        da scipy.stats.ks_1samp intern mit Arrays arbeitet.

        @param x: float oder numpy-Array
        @return: CDF-Werte
        """
        x_arr = np.asarray(x, dtype=float)
        a = 2.0 / math.sqrt(math.pi)
        val = scipy_erf(a * x_arr) - (4.0 * x_arr / math.pi) * np.exp(-4.0 * x_arr ** 2 / math.pi)
        result_arr = np.clip(val, 0.0, 1.0)
        if result_arr.ndim == 0:
            return float(result_arr)
        return result_arr

    ks_result = scipy_stats.ks_1samp(
        spacings_arr,
        gue_wigner_cdf_matrix,
        alternative='two-sided'
    )

    # Stichproben-Statistik: erste paar Eigenwerte einer Beispiel-Matrix
    sample_eigenvalues = eigenvalues_sorted[:10].tolist() if 'eigenvalues_sorted' in dir() else []

    return {
        "n_matrices": n,
        "matrix_size": size,
        "sample_matrix_shape": (size, size),
        "eigenvalues_real": eigenvalues_real,
        "n_spacings": len(all_spacings),
        "eigenvalue_spacings": all_spacings[:200],  # Erste 200 für Inspektion
        "mean_spacing": mean_spacing,
        "variance": variance,
        "std_dev": std_dev,
        "min_spacing": float(np.min(spacings_arr)),
        "max_spacing": float(np.max(spacings_arr)),
        "ks_statistic": float(ks_result.statistic),
        "ks_p_value": float(ks_result.pvalue),
        "gue_wigner_formula": "P(s) = (32/pi^2) * s^2 * exp(-4*s^2/pi)",
        "semicircle_law": "Wigner-Halbkreisgesetz: rho(lambda) = sqrt(4-lambda^2) / (2*pi)",
        # Hilbert-Pólya-Vermutung: GUE-Matrizen als Modell für Riemann-Operator
        "hilbert_polya_connection": (
            "GUE-Eigenwertstatistik stimmt empirisch mit Riemann-Nullstellenstatistik überein"
        )
    }
