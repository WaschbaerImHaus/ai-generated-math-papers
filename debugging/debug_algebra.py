#!/usr/bin/env python3
"""
@file debug_algebra.py
@brief Debug-Skript für das Algebra-Modul.
@description
    Führt exemplarische Berechnungen mit dem Algebra-Modul durch
    und gibt Zwischenschritte und Ergebnisse detailliert aus.
    Dient zum schnellen Überprüfen von Algebra-Funktionen ohne Tests laufen zu lassen.

@author Kurt Ingwer
@date 2026-03-07
@lastModified 2026-03-07
"""

import sys
import os

# Projektverzeichnis zum Suchpfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from algebra import (
    Polynomial, solve_linear, solve_quadratic,
    gcd, lcm, extended_gcd, mod_inverse,
    is_prime, prime_factorization, euler_phi
)

def debug_polynomial():
    """Polynome debuggen."""
    print("=" * 50)
    print("POLYNOME")
    print("=" * 50)

    # p(x) = 3x^2 + 2x - 1
    p = Polynomial([3, 2, -1])
    print(f"p(x) = {p}")
    print(f"p(0) = {p.evaluate(0)}")
    print(f"p(1) = {p.evaluate(1)}")
    print(f"p(-1) = {p.evaluate(-1)}")
    print(f"p'(x) = {p.derivative()}")

    # q(x) = x + 1
    q = Polynomial([1, 1])
    print(f"\nq(x) = {q}")
    print(f"p + q = {p + q}")
    print(f"p * q = {p * q}")


def debug_equation_solver():
    """Gleichungslöser debuggen."""
    print("\n" + "=" * 50)
    print("GLEICHUNGSLÖSER")
    print("=" * 50)

    # Lineare Gleichung: 2x + 4 = 0  →  x = -2
    sol = solve_linear(2, 4)
    print(f"2x + 4 = 0  →  x = {sol}")

    # Quadratische Gleichung: x^2 - 5x + 6 = 0  →  x=2, x=3
    roots = solve_quadratic(1, -5, 6)
    print(f"x² - 5x + 6 = 0  →  {roots}")

    # Komplexe Wurzeln: x^2 + 1 = 0
    roots_complex = solve_quadratic(1, 0, 1)
    print(f"x² + 1 = 0  →  {roots_complex}")


def debug_number_theory():
    """Zahlentheorie debuggen."""
    print("\n" + "=" * 50)
    print("ZAHLENTHEORIE")
    print("=" * 50)

    print(f"gcd(48, 18) = {gcd(48, 18)}  (erwartet: 6)")
    print(f"lcm(4, 6) = {lcm(4, 6)}  (erwartet: 12)")

    g, x, y = extended_gcd(35, 15)
    print(f"extended_gcd(35, 15): g={g}, Bezout: 35*{x} + 15*{y} = {35*x + 15*y}")

    inv = mod_inverse(3, 7)
    print(f"mod_inverse(3, 7) = {inv}  (Check: 3*{inv} mod 7 = {(3*inv) % 7})")

    print(f"is_prime(17) = {is_prime(17)}")
    print(f"is_prime(18) = {is_prime(18)}")
    print(f"prime_factorization(360) = {prime_factorization(360)}")
    print(f"euler_phi(12) = {euler_phi(12)}  (erwartet: 4)")


if __name__ == '__main__':
    debug_polynomial()
    debug_equation_solver()
    debug_number_theory()
    print("\n✓ Debug-Lauf abgeschlossen")
