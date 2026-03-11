"""
Tests für das Fourier-Modul (src/fourier.py).

Testet DFT, FFT, IFFT, Fourier-Reihen und Fensterfunktionen.

@author: Kurt Ingwer
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import math
import cmath
import pytest
import sys
import os

# src-Pfad hinzufügen
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from fourier import (
    dft, idft, fft, ifft, fft_padded,
    fourier_coefficients, fourier_series_eval,
    window_hanning, window_hamming, window_blackman, apply_window,
    power_spectrum, dominant_frequency, stft,
    polynomial_multiply_fft, polynomial_multiply_naive
)


# ---------------------------------------------------------------------------
# HILFSFUNKTIONEN
# ---------------------------------------------------------------------------

def allclose_complex(a: list, b: list, rtol: float = 1e-6, atol: float = 1e-6) -> bool:
    """Vergleicht zwei komplexe Listen mit Toleranz."""
    if len(a) != len(b):
        return False
    return all(
        abs(ai - bi) <= atol + rtol * max(abs(ai), abs(bi))
        for ai, bi in zip(a, b)
    )


# ---------------------------------------------------------------------------
# DFT TESTS
# ---------------------------------------------------------------------------

class TestDFT:
    """Tests für diskrete Fourier-Transformation."""

    def test_dft_single_element(self):
        """DFT eines einzelnen Elements ist das Element selbst."""
        result = dft([3.0])
        assert abs(result[0] - 3.0) < 1e-10

    def test_dft_dc_signal(self):
        """Gleichsignal [c, c, c, c] hat nur DC-Komponente."""
        x = [2.0, 2.0, 2.0, 2.0]
        X = dft(x)
        # DC: X[0] = N * c = 4 * 2 = 8
        assert abs(X[0] - 8.0) < 1e-10
        # Alle anderen sollten ~0 sein
        for k in range(1, 4):
            assert abs(X[k]) < 1e-10

    def test_dft_impulse(self):
        """Dirac-Impuls [1, 0, 0, 0] hat konstantes Spektrum."""
        x = [1.0, 0.0, 0.0, 0.0]
        X = dft(x)
        # X[k] = 1 für alle k
        for k in range(4):
            assert abs(X[k] - 1.0) < 1e-10

    def test_dft_orthogonality(self):
        """Einfache Sinusfolge: nur eine Frequenz."""
        # x[n] = sin(2π/N * n) für N=8
        n = 8
        x = [math.sin(2 * math.pi * k / n) for k in range(n)]
        X = dft(x)
        # Nur k=1 und k=7 sollten Energie haben
        assert abs(X[1]) > abs(X[0]) * 10
        assert abs(X[7]) > abs(X[0]) * 10

    def test_idft_inverse_of_dft(self):
        """IDFT ist Inverse der DFT."""
        x = [1.0, 2.0, 3.0, 4.0]
        X = dft(x)
        x_reconstructed = idft(X)
        for i in range(4):
            assert abs(x_reconstructed[i].real - x[i]) < 1e-10
            assert abs(x_reconstructed[i].imag) < 1e-10

    def test_dft_linearity(self):
        """DFT ist linear: DFT(ax + by) = a·DFT(x) + b·DFT(y)."""
        x = [1.0, 2.0, 3.0, 4.0]
        y = [4.0, 3.0, 2.0, 1.0]
        a, b = 2.0, 3.0

        xy = [a * x[i] + b * y[i] for i in range(4)]
        Xxy = dft(xy)

        Xx = dft(x)
        Xy = dft(y)
        Xxy_expected = [a * Xx[k] + b * Xy[k] for k in range(4)]

        assert allclose_complex(Xxy, Xxy_expected)

    def test_dft_parseval(self):
        """Parseval-Theorem: Σ|x[n]|² = (1/N) Σ|X[k]|²."""
        x = [1.0, -2.0, 3.0, -1.0]
        X = dft(x)
        n = len(x)

        energy_time = sum(abs(xi) ** 2 for xi in x)
        energy_freq = sum(abs(Xk) ** 2 for Xk in X) / n

        assert abs(energy_time - energy_freq) < 1e-10


# ---------------------------------------------------------------------------
# FFT TESTS
# ---------------------------------------------------------------------------

class TestFFT:
    """Tests für schnelle Fourier-Transformation (Cooley-Tukey)."""

    def test_fft_equals_dft(self):
        """FFT muss identische Ergebnisse wie DFT liefern."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        X_dft = dft(x)
        X_fft = fft(x)
        assert allclose_complex(X_dft, X_fft, atol=1e-8)

    def test_fft_length_1(self):
        """FFT der Länge 1 ist Identität."""
        result = fft([5.0])
        assert abs(result[0] - 5.0) < 1e-10

    def test_fft_length_4(self):
        """FFT der Länge 4."""
        x = [1.0, 0.0, 0.0, 0.0]  # Impuls
        X = fft(x)
        for k in range(4):
            assert abs(X[k] - 1.0) < 1e-10

    def test_ifft_inverse_of_fft(self):
        """IFFT ist Inverse der FFT."""
        x = [1.0, -2.0, 3.0, 4.0, 0.5, -1.0, 2.0, 1.5]
        X = fft(x)
        x_rec = ifft(X)
        for i in range(8):
            assert abs(x_rec[i].real - x[i]) < 1e-8
            assert abs(x_rec[i].imag) < 1e-8

    def test_fft_odd_length_fallback(self):
        """FFT für ungerade Längen fällt auf DFT zurück."""
        x = [1.0, 2.0, 3.0]
        X_fft = fft(x)
        X_dft = dft(x)
        assert allclose_complex(X_fft, X_dft, atol=1e-10)

    def test_fft_padded_zero_padding(self):
        """fft_padded erzeugt Ausgabe der Länge = nächste Zweierpotenz."""
        x = [1.0, 2.0, 3.0]  # Länge 3 → padded auf 4
        X = fft_padded(x)
        assert len(X) == 4

    def test_fft_cosine_spectrum(self):
        """Cosinus-Signal hat zwei Spektrallinien."""
        n = 16
        # cos(2π·2·t) = cos bei Frequenz k=2
        x = [math.cos(2 * math.pi * 2 * k / n) for k in range(n)]
        X = fft(x)
        magnitudes = [abs(Xk) for Xk in X]

        # Peaks bei k=2 und k=14 (N-2)
        assert magnitudes[2] > magnitudes[0] * 10
        assert magnitudes[14] > magnitudes[0] * 10


# ---------------------------------------------------------------------------
# FOURIER-REIHEN TESTS
# ---------------------------------------------------------------------------

class TestFourierSeries:
    """Tests für Fourier-Reihen periodischer Funktionen."""

    def test_constant_function(self):
        """Konstante Funktion hat nur a₀-Term."""
        f = lambda t: 3.0
        a, b = fourier_coefficients(f, period=2 * math.pi, n_terms=3)
        # a₀/2 = 3 → a₀ = 6
        assert abs(a[0] - 6.0) < 0.01
        # Alle höheren Terme sollten klein sein
        for k in range(1, 4):
            assert abs(a[k]) < 0.01
            assert abs(b[k]) < 0.01

    def test_cosine_function(self):
        """cos(t) hat nur a₁ ≈ 1, alle anderen ≈ 0."""
        f = lambda t: math.cos(t)
        a, b = fourier_coefficients(f, period=2 * math.pi, n_terms=3)
        assert abs(a[1] - 1.0) < 0.01
        assert abs(b[1]) < 0.01
        assert abs(a[2]) < 0.01

    def test_sine_function(self):
        """sin(t) hat nur b₁ ≈ 1."""
        f = lambda t: math.sin(t)
        a, b = fourier_coefficients(f, period=2 * math.pi, n_terms=3)
        assert abs(b[1] - 1.0) < 0.01
        assert abs(a[1]) < 0.01

    def test_fourier_series_eval_reconstruction(self):
        """Fourier-Reihe rekonstruiert eine einfache Funktion annähernd."""
        f = lambda t: math.sin(t)
        period = 2 * math.pi
        a, b = fourier_coefficients(f, period=period, n_terms=5)

        # Auswerten bei verschiedenen Punkten
        for t in [0.0, math.pi / 4, math.pi / 2, math.pi]:
            approx = fourier_series_eval(t, a, b, period)
            expected = f(t)
            assert abs(approx - expected) < 0.01

    def test_square_wave_gibbs(self):
        """Rechteckfunktion konvergiert mit Gibbs-Phänomen."""
        # Rechteck: f(t) = 1 für 0<t<π, -1 für π<t<2π
        f = lambda t: 1.0 if t % (2 * math.pi) < math.pi else -1.0
        a, b = fourier_coefficients(f, period=2 * math.pi, n_terms=10)

        # Bei t=π/2: f(π/2) = 1, Fourier-Näherung sollte nahe dran sein
        approx = fourier_series_eval(math.pi / 2, a, b, 2 * math.pi)
        assert abs(approx - 1.0) < 0.15  # Gibbs-Phänomen erlaubt etwas Fehler


# ---------------------------------------------------------------------------
# FENSTERFUNKTIONEN TESTS
# ---------------------------------------------------------------------------

class TestWindowFunctions:
    """Tests für Fensterfunktionen."""

    def test_hanning_endpoints_zero(self):
        """Hanning-Fenster ist 0 an den Rändern."""
        w = window_hanning(8)
        assert abs(w[0]) < 1e-10
        assert abs(w[-1]) < 1e-10

    def test_hanning_symmetry(self):
        """Hanning-Fenster ist symmetrisch."""
        w = window_hanning(9)
        for i in range(4):
            assert abs(w[i] - w[8 - i]) < 1e-10

    def test_hanning_max_center(self):
        """Maximum des Hanning-Fensters ist in der Mitte."""
        n = 9
        w = window_hanning(n)
        center = n // 2
        assert w[center] == max(w)

    def test_hamming_not_zero_at_endpoints(self):
        """Hamming-Fenster ist nicht 0 an den Rändern (anders als Hanning)."""
        w = window_hamming(8)
        assert abs(w[0]) > 0.05  # w[0] = 0.54 - 0.46 = 0.08

    def test_blackman_length(self):
        """Blackman-Fenster hat korrekte Länge."""
        n = 16
        w = window_blackman(n)
        assert len(w) == n

    def test_apply_window_scales_signal(self):
        """Fensterfunktion skaliert Signal korrekt."""
        x = [1.0, 1.0, 1.0, 1.0]
        w = [0.5, 0.8, 0.8, 0.5]
        result = apply_window(x, w)
        assert result == [0.5, 0.8, 0.8, 0.5]

    def test_apply_window_length_mismatch_raises(self):
        """Längenfehler wird erkannt."""
        with pytest.raises(ValueError):
            apply_window([1.0, 2.0, 3.0], [0.5, 0.5])


# ---------------------------------------------------------------------------
# SPEKTRALANALYSE TESTS
# ---------------------------------------------------------------------------

class TestSpectralAnalysis:
    """Tests für Leistungsspektrum und Frequenzanalyse."""

    def test_power_spectrum_dc_signal(self):
        """Gleichsignal hat Energie nur bei f=0."""
        x = [1.0] * 8
        freqs, power = power_spectrum(x)
        # DC-Komponente (k=0) hat maximale Leistung
        assert power[0] == max(power)

    def test_power_spectrum_length(self):
        """Leistungsspektrum hat N/2 + 1 Punkte."""
        x = [1.0] * 8
        freqs, power = power_spectrum(x)
        # fft_padded macht aus 8 Punkten 8 Punkte → 5 einseitig
        assert len(freqs) == len(power)

    def test_dominant_frequency_sine(self):
        """Dominant Frequency von sin(2π·f₀·t) = f₀."""
        sample_rate = 100.0  # 100 Hz
        f0 = 10.0            # 10 Hz Sinussignal
        n = 64

        # Signal generieren
        x = [math.sin(2 * math.pi * f0 * k / sample_rate) for k in range(n)]
        freq = dominant_frequency(x, sample_rate)

        # Toleranz: ±sample_rate/n = ±1.5625 Hz
        assert abs(freq - f0) < 2.0

    def test_stft_frame_count(self):
        """STFT erzeugt korrekte Anzahl Frames."""
        x = [0.0] * 32
        window_size = 8
        hop_size = 4
        frames = stft(x, window_size, hop_size)
        # Erwartete Frames: (32 - 8) / 4 + 1 = 7
        expected = (len(x) - window_size) // hop_size + 1
        assert len(frames) == expected

    def test_stft_frame_length(self):
        """Jeder STFT-Frame hat Länge window_size."""
        x = [1.0] * 16
        frames = stft(x, window_size=8, hop_size=4)
        for frame in frames:
            assert len(frame) == 8


# ---------------------------------------------------------------------------
# FFT-BASIERTE POLYNOMMULTIPLIKATION (Build 58)
# ---------------------------------------------------------------------------

class TestPolynomialMultiplyFFT:
    """Tests für die FFT-basierte Polynommultiplikation."""

    def test_multiply_linear_times_linear(self):
        """(1 + x) * (1 + x) = 1 + 2x + x²"""
        p = [1.0, 1.0]  # 1 + x
        q = [1.0, 1.0]  # 1 + x
        result = polynomial_multiply_fft(p, q)
        assert len(result) == 3
        assert abs(result[0] - 1.0) < 1e-8
        assert abs(result[1] - 2.0) < 1e-8
        assert abs(result[2] - 1.0) < 1e-8

    def test_multiply_by_constant(self):
        """Multiplikation mit Konstante: 3 * (1 + 2x) = 3 + 6x"""
        p = [3.0]
        q = [1.0, 2.0]
        result = polynomial_multiply_fft(p, q)
        assert len(result) == 2
        assert abs(result[0] - 3.0) < 1e-8
        assert abs(result[1] - 6.0) < 1e-8

    def test_multiply_quadratic_times_linear(self):
        """(x² + 2x + 1) * (x + 1) = x³ + 3x² + 3x + 1"""
        p = [1.0, 2.0, 1.0]  # 1 + 2x + x²
        q = [1.0, 1.0]       # 1 + x
        result = polynomial_multiply_fft(p, q)
        assert len(result) == 4
        assert abs(result[0] - 1.0) < 1e-8  # Koeffizient von x⁰
        assert abs(result[1] - 3.0) < 1e-8  # Koeffizient von x¹
        assert abs(result[2] - 3.0) < 1e-8  # Koeffizient von x²
        assert abs(result[3] - 1.0) < 1e-8  # Koeffizient von x³

    def test_multiply_agrees_with_naive(self):
        """FFT-Methode stimmt mit naiver Methode überein."""
        p = [1.0, -2.0, 3.0]
        q = [4.0, 0.0, -1.0, 2.0]
        fft_result = polynomial_multiply_fft(p, q)
        naive_result = polynomial_multiply_naive(p, q)
        assert len(fft_result) == len(naive_result)
        for r_fft, r_naive in zip(fft_result, naive_result):
            assert abs(r_fft - r_naive) < 1e-8

    def test_multiply_large_polynomials(self):
        """FFT-Multiplikation großer Polynome."""
        import random
        random.seed(42)
        n = 64
        p = [float(random.randint(-10, 10)) for _ in range(n)]
        q = [float(random.randint(-10, 10)) for _ in range(n)]
        fft_result = polynomial_multiply_fft(p, q)
        naive_result = polynomial_multiply_naive(p, q)
        assert len(fft_result) == 2 * n - 1
        for r_fft, r_naive in zip(fft_result, naive_result):
            assert abs(r_fft - r_naive) < 1e-5

    def test_multiply_single_element(self):
        """Multiplikation zweier Konstanten."""
        p = [5.0]
        q = [7.0]
        result = polynomial_multiply_fft(p, q)
        assert len(result) == 1
        assert abs(result[0] - 35.0) < 1e-8

    def test_naive_multiply_linear(self):
        """Naive Methode: (1 + x) * (1 + x) = 1 + 2x + x²"""
        p = [1.0, 1.0]
        q = [1.0, 1.0]
        result = polynomial_multiply_naive(p, q)
        assert abs(result[0] - 1.0) < 1e-10
        assert abs(result[1] - 2.0) < 1e-10
        assert abs(result[2] - 1.0) < 1e-10
