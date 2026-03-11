"""
@file modular_forms_theta.py
@brief Theta-Reihen, Jacobi-Theta-Funktionen und Dedekind-Eta für Modulformen.
@description
    Implementiert Theta-Reihen und verwandte Funktionen:
    - Jacobi-Theta-Funktion ϑ_3(0|τ)
    - Theta-Transformationsformel (Jacobiische Imagination)
    - Darstellungsanzahlen r_k(n) als Summe von k Quadraten (via Theta-Reihen)
    - Jacobi-Dreifachprodukt (Jacobi 1829)
    - Dedekind-Eta-Funktion η(τ)

    Theta-Reihen sind fundamentale Modulformen halben Gewichts (Gewicht 1/2).
    Sie verbinden:
    - Analytische Zahlentheorie (Darstellungsanzahlen)
    - Physik (Quantenstring-Theorie, Wärmeleitungsgleichung)
    - Elliptische Kurven (Verbindung via Modulformen)

    Die Theta-Transformationsformel:
        θ(-1/τ) = sqrt(-iτ) · θ(τ)
    ist äquivalent zur Poissonschen Summenformel und war das Werkzeug, das
    Riemann für die Funktionalgleichung der Riemann-Zeta-Funktion verwendete.

    Die Dedekind-Eta-Funktion η(τ) erfüllt:
        Δ(τ) = η(τ)^{24}
    (Verbindung zur Diskriminantenform)

@author Michael Fuhrmann
@version 1.0
@since 2026-03-11
@lastModified 2026-03-11
"""

import math
import cmath


# ===========================================================================
# THETA-REIHEN
# ===========================================================================

def theta_function(z: complex, n_max: int = 50) -> complex:
    """
    Berechnet die Jacobi-Theta-Funktion ϑ_3(0|τ).

    Definition (Jacobi-Theta-Funktion dritten Typs bei z=0):
        θ(τ) = ϑ_3(0|τ) = Σ_{n=-∞}^{∞} q^{n²}   mit q = e^{πiτ}

    Äquivalent: θ(τ) = 1 + 2·Σ_{n=1}^∞ q^{n²} (da n und -n denselben Beitrag liefern)

    Konvergenz: Für Im(τ) > 0 gilt |q| = e^{-π·Im(τ)} < 1,
    daher konvergiert die Reihe exponentiell schnell.

    Anwendungen:
    - θ(τ)^k = Σ_{n=0}^∞ r_k(n) q^n liefert Darstellungsanzahlen r_k(n)
    - θ ist eine Modulform halben Gewichts (Gewicht 1/2)
    - Verbindung zu elliptischen Funktionen und Modulformen

    Transformationsformel (Jacobiische Imagination):
        θ(-1/τ) = sqrt(-iτ) · θ(τ)

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_max: Maximaler Summationsindex (symmetrisch: -n_max bis n_max)
    @return: Wert der Theta-Funktion als komplexe Zahl
    @raises ValueError: Wenn Im(z) ≤ 0
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if z.imag <= 0:
        raise ValueError(f"τ muss in der oberen Halbebene liegen (Im(τ) > 0), erhalten: {z}")

    # q = e^{πiτ}  (beachte: π, nicht 2π !)
    q = cmath.exp(1j * math.pi * z)

    # θ(τ) = Σ_{n=-N}^{N} q^{n²} = 1 + 2·Σ_{n=1}^{N} q^{n²}
    # Nutze symmetrische Form für Effizienz (n und -n tragen gleich bei)
    total = complex(1.0)  # n=0 Beitrag
    for n in range(1, n_max + 1):
        # q^{n²} berechnen
        term = q ** (n * n)
        # Konvergenzcheck: bei sehr kleinen Werten abbrechen
        if abs(term) < 1e-300:
            break
        # Faktor 2 wegen n und -n (gleicher Beitrag da n²=(-n)²)
        total += 2.0 * term

    return total


def theta_transformation(z: complex) -> dict[str, object]:
    """
    Verifiziert die Jacobi-Theta-Transformationsformel numerisch.

    Die Transformationsformel lautet:
        θ(-1/τ) = sqrt(-iτ) · θ(τ)

    Dies ist eine fundamentale Symmetrie der Theta-Funktion unter der
    S-Transformation (τ ↦ -1/τ) der modularen Gruppe.

    Herleitung via Poissonsche Summenformel:
        Σ_{n} e^{-πn²/t} = √t · Σ_{n} e^{-πn²t}
    (Riemann verwendete diese Formel für seine zeta-Funktion!)

    @param z: Punkt τ in der oberen Halbebene
    @return: Dictionary mit:
             'theta_z': θ(τ)
             'theta_minus_inv': θ(-1/τ)
             'transformation_ratio': θ(-1/τ) / θ(τ)
             'expected': sqrt(-iτ) (erwartetes Verhältnis)
             'verified': True wenn |ratio - expected| / |expected| < 1e-4
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Berechne θ(τ)
    theta_z = theta_function(z)

    # Berechne θ(-1/τ) – der transformierte Argument ist -1/τ
    z_transformed = -1.0 / z
    theta_minus_inv = theta_function(z_transformed)

    # Erwartetes Verhältnis: sqrt(-iτ)
    # sqrt(-iτ) = sqrt(-i) · sqrt(τ) mit Hauptzweig
    expected = cmath.sqrt(-1j * z)

    # Tatsächliches Verhältnis
    if abs(theta_z) < 1e-300:
        # Sonderfall: θ(τ) ≈ 0 (sollte nicht vorkommen für Im(τ) > 0)
        transformation_ratio = complex(float('nan'))
        verified = False
    else:
        transformation_ratio = theta_minus_inv / theta_z
        # Relativer Fehler
        rel_error = abs(transformation_ratio - expected) / max(abs(expected), 1e-15)
        verified = rel_error < 1e-4

    return {
        'theta_z': theta_z,
        'theta_minus_inv': theta_minus_inv,
        'transformation_ratio': transformation_ratio,
        'expected': expected,
        'verified': verified
    }


def sum_of_squares_theta(n: int, k: int = 2) -> int:
    """
    Berechnet r_k(n): Anzahl der Darstellungen von n als Summe von k Quadraten.

    Definition:
        r_k(n) = #{(x_1,...,x_k) ∈ Z^k : x_1² + ... + x_k² = n}
    (negative Zahlen und verschiedene Reihenfolgen zählen separat)

    Theoretisch via Theta-Reihen:
        θ(τ)^k = (Σ_{n} q^{n²})^k = Σ_{m=0}^∞ r_k(m) q^m
    d.h. r_k(m) ist der m-te Koeffizient von θ^k.

    Für k=2 gilt die Formel von Jacobi (1829):
        r_2(n) = 4 · Σ_{d|n} χ(d)
    wobei χ(d) = 0 für gerades d, χ(d) = (-1)^{(d-1)/2} für ungerades d.
    (Dirichlet-Charakter mod 4)

    Bekannte Werte:
        r_2(0) = 1, r_2(1) = 4, r_2(2) = 4, r_2(4) = 4
        r_2(5) = 8 (5 = 1²+2² = 2²+1² = (-1)²+2² = 2²+(-1)² = 1²+(-2)² = ...)
        r_2(25) = 12

    Diese Implementierung nutzt Brute-Force für alle k (universell und korrekt).

    @param n: Nicht-negative ganze Zahl
    @param k: Anzahl der Quadrate (Standard: 2)
    @return: r_k(n) als nicht-negative ganze Zahl
    @raises ValueError: Wenn n < 0 oder k < 1
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if n < 0:
        raise ValueError(f"n muss ≥ 0 sein, erhalten: {n}")
    if k < 1:
        raise ValueError(f"k muss ≥ 1 sein, erhalten: {k}")

    # Sonderfälle
    if n == 0:
        return 1  # Nur die Nulldarstellung: (0,0,...,0)

    # Maximaler Betrag einer Komponente: |x_i| ≤ floor(sqrt(n))
    max_val = int(math.isqrt(n))

    # Rekursive Zählung via dynamischer Programmierung
    # count[m] = Anzahl Wege, m als Summe von genau j Quadraten darzustellen
    # Iteriere über die k Quadrate
    # Starte mit 1 Weg für Summe=0 (leere Summe)
    current = {0: 1}

    for _ in range(k):
        # Nächste Iteration: addiere ein weiteres Quadrat x² für x ∈ Z
        next_count: dict[int, int] = {}
        for s, cnt in current.items():
            # Iteriere über alle möglichen x (von -max_val bis max_val)
            for x in range(-max_val, max_val + 1):
                new_s = s + x * x
                if new_s <= n:  # Nur Werte ≤ n weiter verfolgen
                    next_count[new_s] = next_count.get(new_s, 0) + cnt
        current = next_count

    return current.get(n, 0)


def jacobi_triple_product(z: complex, q: complex, n_terms: int = 20) -> complex:
    """
    Berechnet das Jacobi-Dreifachprodukt.

    Das Jacobi-Dreifachprodukt (Jacobi, 1829) ist eine der schönsten Identitäten
    der Mathematik:
        Π_{n=1}^∞ (1 - q^{2n})(1 + z·q^{2n-1})(1 + z^{-1}·q^{2n-1})
        = Σ_{n=-∞}^{∞} z^n · q^{n²}

    Voraussetzung: |q| < 1 und z ≠ 0.

    Spezialfälle:
        z = 1: Gibt θ_3(0|τ) = Σ q^{n²} zurück (mit q = e^{πiτ})
        z = -1: Gibt θ_4(0|τ) zurück
        z = q: Gibt θ_2(0|τ) zurück

    Anwendungen:
        - Verbindung zu Eta-Funktion: η(τ)³ = q^{1/8}·Σ_{n} (-1)^n·(2n+1)·q^{n(n+1)/2}
        - Beweis der Pentagonalzahlformel von Euler
        - Partitionentheorie

    @param z: Komplexer Parameter (z ≠ 0)
    @param q: Konvergenzparameter mit |q| < 1
    @param n_terms: Anzahl der Produktterme
    @return: Wert des Dreifachprodukts (Produktseite)
    @raises ValueError: Wenn |q| ≥ 1 oder z = 0
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if abs(q) >= 1.0:
        raise ValueError(f"|q| muss < 1 sein, erhalten: |q| = {abs(q)}")
    if abs(z) < 1e-300:
        raise ValueError("z darf nicht 0 sein")

    # Produktformel: Π_{n=1}^{N} (1-q^{2n})(1+z·q^{2n-1})(1+z^{-1}·q^{2n-1})
    product = complex(1.0)
    z_inv = 1.0 / z

    for n in range(1, n_terms + 1):
        # q^{2n} und q^{2n-1} berechnen
        q2n = q ** (2 * n)
        q2n_m1 = q ** (2 * n - 1)

        # Drei Faktoren des Produkts
        f1 = 1.0 - q2n                   # (1 - q^{2n})
        f2 = 1.0 + z * q2n_m1            # (1 + z·q^{2n-1})
        f3 = 1.0 + z_inv * q2n_m1        # (1 + z^{-1}·q^{2n-1})

        product *= f1 * f2 * f3

        # Konvergenzcheck: wenn q^{2n} sehr klein, abbrechen
        if abs(q2n) < 1e-300:
            break

    return product


def dedekind_eta(z: complex, n_terms: int = 100) -> complex:
    """
    Berechnet die Dedekind-Eta-Funktion η(τ).

    Definition (Dedekind, 1877):
        η(τ) = q^{1/24} · Π_{n=1}^∞ (1 - q^n)   mit q = e^{2πiτ}

    Eigenschaften:
        - η ist eine Modulform halben Gewichts (Gewicht 1/2) mit Multiplikatorsystem
        - Transformationsformel: η(-1/τ) = sqrt(-iτ) · η(τ)
        - Periodizität: η(τ+1) = e^{πi/12} · η(τ)
        - Verbindung zur Delta-Funktion: Δ(τ) = η(τ)^{24}
          (d.h. η(τ) = Δ(τ)^{1/24})
        - Etafunktion taucht in der Physik auf: Quantenstring-Theorie (Partition function)

    Pentagonalzahlformel von Euler:
        Π_{n=1}^∞ (1-q^n) = Σ_{n=-∞}^∞ (-1)^n · q^{n(3n-1)/2}
    (Pentagonalzahlen: 0, 1, 2, 5, 7, 12, 15, ...)

    @param z: Punkt τ in der oberen Halbebene (Im(z) > 0)
    @param n_terms: Anzahl der Produktterme (mehr = genauer)
    @return: Wert η(τ) als komplexe Zahl
    @raises ValueError: Wenn Im(z) ≤ 0
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if z.imag <= 0:
        raise ValueError(f"τ muss in der oberen Halbebene liegen (Im(τ) > 0), erhalten: {z}")

    # q = e^{2πiτ}  (ganzer Kreis!)
    q = cmath.exp(2j * math.pi * z)

    # Vorfaktor: q^{1/24} = e^{2πiτ/24} = e^{πiτ/12}
    q_124 = cmath.exp(2j * math.pi * z / 24.0)

    # Unendliches Produkt: Π_{n=1}^{N} (1 - q^n)
    product = complex(1.0)
    q_power = q  # q^n, beginnt bei q^1

    for n in range(1, n_terms + 1):
        if abs(q_power) < 1e-300:
            break
        product *= (1.0 - q_power)
        q_power *= q  # q^{n+1}

    # η(τ) = q^{1/24} · Π(1-q^n)
    return q_124 * product
