"""
Fourier-Analysis – DFT, FFT, Fourier-Reihen, Fensterfunktionen.

Dieses Modul implementiert die diskrete und schnelle Fourier-Transformation
sowie Fourier-Reihen für periodische Funktionen:
    - DFT (Diskrete Fourier-Transformation, O(n²))
    - FFT (Cooley-Tukey, O(n log n))
    - Inverse FFT (IFFT)
    - Fourier-Reihen (trigonometrische Koeffizienten)
    - Fensterfunktionen (Hanning, Hamming, Blackman, Kaiser)
    - Spektralanalyse und Leistungsdichtespektrum

Mathematischer Hintergrund:
    Fourier-Transformation (kontinuierlich):
        F(ω) = ∫_{-∞}^{∞} f(t) e^{-iωt} dt
    Diskrete Fourier-Transformation:
        X[k] = Σ_{n=0}^{N-1} x[n] · e^{-2πi·kn/N}
    Inverse DFT:
        x[n] = (1/N) · Σ_{k=0}^{N-1} X[k] · e^{2πi·kn/N}

Verbindung zur Zahlentheorie:
    Die Fourier-Transformation ist entscheidend für Riemanns Explizite Formel
    und die Analyse der Zeta-Nullstellen über den Poisson-Summationsformel.

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import math
import cmath
import numpy as np
from typing import Callable, Optional


# ===========================================================================
# DISKRETE FOURIER-TRANSFORMATION (DFT)
# ===========================================================================

def dft(x: list) -> list:
    """
    Berechnet die diskrete Fourier-Transformation (DFT).

    Definition:
        X[k] = Σ_{n=0}^{N-1} x[n] · e^{-2πi·kn/N}

    Komplexität: O(N²) – für kleine N oder Lernzwecke.
    Für große N besser fft() verwenden (O(N log N)).

    @param x: Eingangsfolge (real oder komplex)
    @return: Fourier-Koeffizienten X[0..N-1]
    @lastModified: 2026-03-08
    """
    n = len(x)
    # Twiddle-Faktor: ω = e^{-2πi/N}
    omega = cmath.exp(-2j * cmath.pi / n)

    result = []
    for k in range(n):
        # X[k] = Σ x[n] · ω^(kn)
        xk = sum(x[j] * (omega ** (k * j)) for j in range(n))
        result.append(xk)

    return result


def idft(x: list) -> list:
    """
    Berechnet die inverse diskrete Fourier-Transformation (IDFT).

    Definition:
        x[n] = (1/N) · Σ_{k=0}^{N-1} X[k] · e^{2πi·kn/N}

    Trick: IDFT = konjugierte DFT / N

    @param x: Fourier-Koeffizienten X[0..N-1]
    @return: Rekonstruierte Zeitfolge x[0..N-1]
    @lastModified: 2026-03-08
    """
    n = len(x)
    # Konjugiere, DFT anwenden, konjugiere und normiere → IDFT
    x_conj = [v.conjugate() for v in x]
    dft_conj = dft(x_conj)
    return [v.conjugate() / n for v in dft_conj]


# ===========================================================================
# SCHNELLE FOURIER-TRANSFORMATION (FFT) – COOLEY-TUKEY
# ===========================================================================

def fft(x: list) -> list:
    """
    Schnelle Fourier-Transformation (FFT) nach Cooley-Tukey (radix-2).

    Algorithmus (rekursiv):
        X[k] = X_even[k] + ω^k · X_odd[k]   (k < N/2)
        X[k+N/2] = X_even[k] - ω^k · X_odd[k]

    Komplexität: O(N log N) für N = 2^m (Zweierpotenz).
    Für N ≠ 2^m: Nullauffüllung (Zero-Padding) auf nächste Zweierpotenz.

    Historisch: Cooley & Tukey 1965 (wiederentdeckt, da Gauß ~1805 bereits kannte).

    @param x: Eingangsfolge, Länge sollte Zweierpotenz sein
    @return: Fourier-Koeffizienten X[0..N-1]
    @lastModified: 2026-03-08
    """
    n = len(x)

    # Basisfall: N=1
    if n == 1:
        return [complex(x[0])]

    # Wenn N keine Zweierpotenz: Fallback auf DFT
    if n % 2 != 0:
        return dft(x)

    # Teile-und-Herrsche: Gerade und ungerade Indizes trennen
    x_even = fft(x[::2])   # x[0], x[2], x[4], ...
    x_odd  = fft(x[1::2])  # x[1], x[3], x[5], ...

    # Twiddle-Faktoren ω^k = e^{-2πik/N}
    twiddle = [cmath.exp(-2j * cmath.pi * k / n) for k in range(n // 2)]

    # Butterfly-Kombination
    result = [0] * n
    for k in range(n // 2):
        t = twiddle[k] * x_odd[k]
        result[k]           = x_even[k] + t  # Obere Hälfte
        result[k + n // 2]  = x_even[k] - t  # Untere Hälfte (Vorzeichen dreht)

    return result


def ifft(x: list) -> list:
    """
    Inverse schnelle Fourier-Transformation (IFFT).

    Benutzt den Trick: IFFT(X) = konj(FFT(konj(X))) / N
    Dadurch ist nur eine FFT-Implementierung nötig.

    @param x: Fourier-Koeffizienten X[0..N-1]
    @return: Rekonstruierte Zeitfolge x[0..N-1]
    @lastModified: 2026-03-08
    """
    n = len(x)
    x_conj = [v.conjugate() for v in x]
    fft_conj = fft(x_conj)
    return [v.conjugate() / n for v in fft_conj]


def fft_padded(x: list) -> list:
    """
    FFT mit automatischer Nullauffüllung auf nächste Zweierpotenz (Zero-Padding).

    Zero-Padding verbessert die Frequenzauflösung im Spektrum.
    Die Ausgabe hat die Länge der nächsten Zweierpotenz ≥ len(x).

    @param x: Eingangsfolge (beliebige Länge)
    @return: FFT-Koeffizienten (Länge = nächste Zweierpotenz ≥ N)
    @lastModified: 2026-03-08
    """
    n = len(x)
    # Nächste Zweierpotenz berechnen: 2^⌈log₂(n)⌉
    n_padded = 1
    while n_padded < n:
        n_padded *= 2

    # Mit Nullen auffüllen
    x_padded = list(x) + [0] * (n_padded - n)
    return fft(x_padded)


# ===========================================================================
# FOURIER-REIHEN FÜR PERIODISCHE FUNKTIONEN
# ===========================================================================

def fourier_coefficients(f: Callable, period: float, n_terms: int,
                         n_samples: int = 1000) -> tuple:
    """
    Berechnet Fourier-Koeffizienten einer periodischen Funktion numerisch.

    Fourier-Reihe für T-periodische Funktion f(t):
        f(t) = a₀/2 + Σ_{k=1}^N [aₖ cos(2πkt/T) + bₖ sin(2πkt/T)]

    Koeffizienten via numerische Integration (Trapezregel):
        aₖ = (2/T) ∫_0^T f(t) cos(2πkt/T) dt
        bₖ = (2/T) ∫_0^T f(t) sin(2πkt/T) dt

    @param f: Periodische Funktion
    @param period: Periode T
    @param n_terms: Anzahl der Terme (0..n_terms)
    @param n_samples: Anzahl der Stützpunkte für numerische Integration
    @return: (a_coeffs, b_coeffs) – Listen der cos/sin-Koeffizienten
    @lastModified: 2026-03-08
    """
    # Äquidistante Stützpunkte über eine Periode
    t_vals = [period * k / n_samples for k in range(n_samples)]
    f_vals = [f(t) for t in t_vals]
    dt = period / n_samples

    a_coeffs = []
    b_coeffs = []

    for k in range(n_terms + 1):
        # Numerische Integration via Trapezregel
        omega_k = 2 * math.pi * k / period

        a_k = (2 / period) * sum(
            f_vals[j] * math.cos(omega_k * t_vals[j]) * dt
            for j in range(n_samples)
        )
        b_k = (2 / period) * sum(
            f_vals[j] * math.sin(omega_k * t_vals[j]) * dt
            for j in range(n_samples)
        )

        a_coeffs.append(a_k)
        b_coeffs.append(b_k)

    return a_coeffs, b_coeffs


def fourier_series_eval(t: float, a_coeffs: list, b_coeffs: list,
                        period: float) -> float:
    """
    Wertet eine Fourier-Reihe an einem Punkt t aus.

    f(t) ≈ a₀/2 + Σ_{k=1}^N [aₖ cos(2πkt/T) + bₖ sin(2πkt/T)]

    @param t: Auswertungspunkt
    @param a_coeffs: Kosinus-Koeffizienten [a₀, a₁, ..., aₙ]
    @param b_coeffs: Sinus-Koeffizienten [b₀, b₁, ..., bₙ]
    @param period: Periode T
    @return: Näherungswert der Fourier-Reihe bei t
    @lastModified: 2026-03-08
    """
    n_terms = len(a_coeffs) - 1
    result = a_coeffs[0] / 2  # Gleichanteil

    for k in range(1, n_terms + 1):
        omega_k = 2 * math.pi * k / period
        result += a_coeffs[k] * math.cos(omega_k * t)
        result += b_coeffs[k] * math.sin(omega_k * t)

    return result


# ===========================================================================
# FENSTERFUNKTIONEN (für spektrale Analyse)
# ===========================================================================

def window_hanning(n: int) -> list:
    """
    Hanning-Fenster (raised cosine): w[k] = 0.5(1 - cos(2πk/N)).

    Reduziert Leck-Effekte (spectral leakage) bei nicht-periodischen Signalen.
    Guter Kompromiss zwischen Hauptkeulen-Breite und Seitenkeulendämpfung.

    @param n: Fensterlänge
    @return: Fensterfunktion als Liste
    @lastModified: 2026-03-08
    """
    return [0.5 * (1 - math.cos(2 * math.pi * k / (n - 1))) for k in range(n)]


def window_hamming(n: int) -> list:
    """
    Hamming-Fenster: w[k] = 0.54 - 0.46·cos(2πk/N).

    Etwas bessere Seitenkeulendämpfung als Hanning (42 dB vs 31 dB).

    @param n: Fensterlänge
    @return: Fensterfunktion als Liste
    @lastModified: 2026-03-08
    """
    return [0.54 - 0.46 * math.cos(2 * math.pi * k / (n - 1)) for k in range(n)]


def window_blackman(n: int) -> list:
    """
    Blackman-Fenster: w[k] = 0.42 - 0.5·cos(2πk/N) + 0.08·cos(4πk/N).

    Sehr gute Seitenkeulendämpfung (74 dB), aber breitere Hauptkeule.
    Gut für Signale mit nahe beieinander liegenden Frequenzen.

    @param n: Fensterlänge
    @return: Fensterfunktion als Liste
    @lastModified: 2026-03-08
    """
    return [
        0.42 - 0.5 * math.cos(2 * math.pi * k / (n - 1))
            + 0.08 * math.cos(4 * math.pi * k / (n - 1))
        for k in range(n)
    ]


def apply_window(x: list, window: list) -> list:
    """
    Wendet eine Fensterfunktion auf eine Signalfolge an.

    Vor der FFT multipliziert man das Signal mit dem Fenster,
    um Randeffekte (Gibbs-Phänomen) zu reduzieren.

    @param x: Eingangssignal
    @param window: Fensterfunktion (gleiche Länge wie x)
    @return: Gefenstertes Signal
    @raises ValueError: Wenn Längen nicht übereinstimmen
    @lastModified: 2026-03-08
    """
    if len(x) != len(window):
        raise ValueError(
            f"Signal ({len(x)}) und Fenster ({len(window)}) müssen gleich lang sein"
        )
    return [x[k] * window[k] for k in range(len(x))]


# ===========================================================================
# SPEKTRALANALYSE
# ===========================================================================

def power_spectrum(x: list, sample_rate: float = 1.0) -> tuple:
    """
    Berechnet das Leistungsdichtespektrum eines Signals.

    Leistungsdichtespektrum (PSD):
        P[k] = |X[k]|² / N²

    Wobei X[k] die FFT-Koeffizienten sind.
    Einseitiges Spektrum (nur positive Frequenzen) für reelle Signale.

    @param x: Reelles Eingangssignal
    @param sample_rate: Abtastrate in Hz (Standard: 1.0)
    @return: (freqs, power) – Frequenzachse und Leistungsspektrum
    @lastModified: 2026-03-08
    """
    n = len(x)
    # FFT berechnen
    x_fft = fft_padded(x)
    n_fft = len(x_fft)

    # Leistungsspektrum: |X[k]|² / N²
    power = [abs(x_fft[k]) ** 2 / n ** 2 for k in range(n_fft // 2 + 1)]

    # Frequenzachse: f[k] = k · sample_rate / N
    freqs = [k * sample_rate / n_fft for k in range(n_fft // 2 + 1)]

    return freqs, power


def dominant_frequency(x: list, sample_rate: float = 1.0) -> float:
    """
    Findet die dominante Frequenz in einem Signal.

    Nützlich für Signalanalyse, Schwingungsanalyse, etc.

    @param x: Eingangssignal
    @param sample_rate: Abtastrate in Hz
    @return: Dominante Frequenz in Hz
    @lastModified: 2026-03-08
    """
    freqs, power = power_spectrum(x, sample_rate)

    # Index des Maximums (außer DC-Komponente bei k=0)
    if len(power) <= 1:
        return 0.0

    max_idx = max(range(1, len(power)), key=lambda k: power[k])
    return freqs[max_idx]


# ===========================================================================
# KURZE FOURIER-TRANSFORMATION (STFT)
# ===========================================================================

def stft(x: list, window_size: int, hop_size: int,
         window_type: str = 'hanning') -> list:
    """
    Kurzzeit-Fourier-Transformation (Short-Time Fourier Transform, STFT).

    Analysiert, wie sich das Frequenzspektrum im Zeitverlauf ändert.
    Das Signal wird in überlappende Fenster unterteilt, auf jedes Fenster
    wird eine FFT angewendet.

    @param x: Eingangssignal
    @param window_size: Fenstergröße (in Samples)
    @param hop_size: Schrittweite zwischen Fenstern (in Samples)
    @param window_type: Fenstertyp ('hanning', 'hamming', 'blackman', 'rect')
    @return: Liste von FFT-Koeffizienten (pro Fenster ein Eintrag)
    @lastModified: 2026-03-08
    """
    n = len(x)

    # Fensterfunktion auswählen
    if window_type == 'hanning':
        win = window_hanning(window_size)
    elif window_type == 'hamming':
        win = window_hamming(window_size)
    elif window_type == 'blackman':
        win = window_blackman(window_size)
    else:
        win = [1.0] * window_size  # Rechteckfenster

    frames = []
    pos = 0
    while pos + window_size <= n:
        # Fenster ausschneiden und mit Fensterfunktion gewichten
        segment = x[pos:pos + window_size]
        windowed = apply_window(segment, win)
        # FFT des gefensterten Segments
        frames.append(fft(windowed))
        pos += hop_size

    return frames
