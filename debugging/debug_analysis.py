#!/usr/bin/env python3
"""
@file debug_analysis.py
@brief Debug-Skript für das Analysis-Modul.
@description
    Führt exemplarische Berechnungen mit dem Analysis-Modul durch
    und zeigt numerische Genauigkeit, Konvergenzverhalten und
    Zwischenschritte detailliert an.

@author Kurt Ingwer
@date 2026-03-07
@lastModified 2026-03-07
"""

import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from analysis import (
    numerical_derivative,
    numerical_integral, newton_raphson, bisection, taylor_series
)


def debug_derivatives():
    """Ableitungen debuggen mit Fehleranalyse."""
    print("=" * 50)
    print("ABLEITUNGEN")
    print("=" * 50)

    # f(x) = sin(x), f'(x) = cos(x)
    f = math.sin
    x = math.pi / 4

    approx = numerical_derivative(f, x)
    exact = math.cos(x)
    error = abs(approx - exact)
    print(f"f(x) = sin(x), x = π/4")
    print(f"  f'(x) numerisch: {approx:.10f}")
    print(f"  f'(x) exakt:     {exact:.10f}")
    print(f"  Fehler:          {error:.2e}")

    # f''(x) = -sin(x)
    approx2 = numerical_derivative(f, x, order=2)
    exact2 = -math.sin(x)
    error2 = abs(approx2 - exact2)
    print(f"\n  f''(x) numerisch: {approx2:.10f}")
    print(f"  f''(x) exakt:     {exact2:.10f}")
    print(f"  Fehler:           {error2:.2e}")


def debug_integration():
    """Numerische Integration debuggen."""
    print("\n" + "=" * 50)
    print("INTEGRATION")
    print("=" * 50)

    # ∫₀^π sin(x) dx = 2
    result = numerical_integral(math.sin, 0, math.pi)
    print(f"∫₀^π sin(x) dx = {result:.10f}  (exakt: 2.0, Fehler: {abs(result-2):.2e})")

    # ∫₀^1 x² dx = 1/3
    result2 = numerical_integral(lambda x: x**2, 0, 1)
    print(f"∫₀^1 x² dx     = {result2:.10f}  (exakt: 0.333..., Fehler: {abs(result2 - 1/3):.2e})")


def debug_root_finding():
    """Nullstellensuche debuggen."""
    print("\n" + "=" * 50)
    print("NULLSTELLENSUCHE")
    print("=" * 50)

    # Newton-Raphson: sqrt(2) als Nullstelle von x^2 - 2
    # Ableitung wird intern numerisch berechnet, kein df nötig
    fn = lambda x: x**2 - 2
    root_nr = newton_raphson(fn, 1.0)
    print(f"Newton-Raphson: √2 ≈ {root_nr:.10f}  (exakt: {math.sqrt(2):.10f})")

    # Bisektion: sin(x) in [3, 4]
    root_bis = bisection(math.sin, 3.0, 4.0)
    print(f"Bisektion: sin(x)=0 in [3,4]: x ≈ {root_bis:.10f}  (exakt: π ≈ {math.pi:.10f})")


if __name__ == '__main__':
    debug_derivatives()
    debug_integration()
    debug_root_finding()
    print("\n✓ Debug-Lauf abgeschlossen")
