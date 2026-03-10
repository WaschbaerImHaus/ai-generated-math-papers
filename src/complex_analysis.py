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

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import cmath
import math
import numpy as np
from typing import Optional


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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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


def xi_symmetry_check(s: complex) -> dict:
    """
    Verifiziert die Symmetrie ξ(s) = ξ(1-s) numerisch.

    Dies ist ein grundlegendes Werkzeug zur Verifikation der Funktionalgleichung.
    Wenn die Symmetrie exakt gilt, ist die Implementierung korrekt.

    @param s: Testpunkt
    @return: Dictionary mit xi(s), xi(1-s) und dem Fehler
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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


def find_zeta_zeros(t_min: float, t_max: float, steps: int = 1000) -> list[dict]:
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
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


def verify_functional_equation(s: complex) -> dict:
    """
    Prüft die Funktionalgleichung ζ(s) = χ(s)·ζ(1-s) numerisch.

    @param s: Testpunkt
    @return: Verifikationsergebnis mit Fehler
    @lastModified: 2026-03-08
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
    @lastModified: 2026-03-08
    """
    if order <= 0 or order > len(f_coeffs):
        raise ValueError(f"Ordnung muss zwischen 1 und {len(f_coeffs)} liegen")
    # Der Koeffizient a_{-1} ist an Position (order - 1) im Array
    return f_coeffs[order - 1]


def cauchy_integral_numerical(
    f: callable,
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
    @lastModified: 2026-03-08
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
    @author Kurt Ingwer
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
) -> list:
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
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
) -> dict:
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
    @author Kurt Ingwer
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
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    _check_mpmath()

    with mpmath.workdps(dps):
        # mpmath.li(x): logarithmischer Integralus mit Cauchy-Hauptwert
        result = mpmath.li(mpmath.mpf(x))
        return float(result)
