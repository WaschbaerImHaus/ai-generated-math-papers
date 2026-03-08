# fourier.py – Fourier-Analysis-Modul

**Autor:** Kurt Ingwer
**Letzte Änderung:** 2026-03-08
**Datei:** `src/fourier.py`

---

## Überblick

Dieses Modul implementiert die Fourier-Analyse – eines der mächtigsten Werkzeuge der angewandten Mathematik. Es umfasst:

| Komponente | Beschreibung |
|-----------|-------------|
| `dft` / `idft` | Diskrete Fourier-Transformation O(N²) |
| `fft` / `ifft` | Schnelle Fourier-Transformation O(N log N) |
| `fft_padded` | FFT mit automatischem Zero-Padding |
| `fourier_coefficients` | Fourier-Koeffizienten periodischer Funktionen |
| `fourier_series_eval` | Auswertung einer Fourier-Reihe |
| `window_hanning/hamming/blackman` | Fensterfunktionen |
| `apply_window` | Anwendung einer Fensterfunktion |
| `power_spectrum` | Leistungsdichtespektrum |
| `dominant_frequency` | Dominante Frequenzkomponente |
| `stft` | Kurzzeit-Fourier-Transformation |

---

## Mathematischer Hintergrund

### Kontinuierliche Fourier-Transformation
```
F(ω) = ∫_{-∞}^{∞} f(t) · e^{-iωt} dt
```

### Diskrete Fourier-Transformation (DFT)
```
X[k] = Σ_{n=0}^{N-1} x[n] · e^{-2πi·kn/N}
```

### Inverse DFT (IDFT)
```
x[n] = (1/N) · Σ_{k=0}^{N-1} X[k] · e^{2πi·kn/N}
```

### Verbindung zur Zahlentheorie
Die Fourier-Transformation ist entscheidend für **Riemanns Explizite Formel** und die Analyse der Zeta-Nullstellen via Poisson-Summationsformel.

---

## DFT und IDFT

### `dft(x) → list`

**Diskrete Fourier-Transformation** (naive Implementierung):

```python
X[k] = Σ_{n=0}^{N-1} x[n] · ω^(kn)    mit ω = e^{-2πi/N}
```

**Komplexität:** O(N²) – geeignet für kleine N oder Lernzwecke.
Für große N: **`fft()`** verwenden (O(N log N)).

### `idft(x) → list`

**Inverse DFT** via Konjugationstrick:
```
IDFT(X) = konj(DFT(konj(X))) / N
```

Dadurch ist nur eine DFT-Implementierung nötig.

---

## FFT (Cooley-Tukey-Algorithmus)

### `fft(x) → list`

**Schnelle Fourier-Transformation** nach Cooley-Tukey (1965), Radix-2-Algorithmus.

### Algorithmus (rekursiv, Teile-und-Herrsche)

```
X_even[k] = FFT(x[0], x[2], x[4], ...)   # Gerade Indizes
X_odd[k]  = FFT(x[1], x[3], x[5], ...)   # Ungerade Indizes

X[k]       = X_even[k] + ω^k · X_odd[k]   (für k < N/2)
X[k + N/2] = X_even[k] - ω^k · X_odd[k]
```

Die **Butterfly-Operation** `(X_even[k] ± ω^k · X_odd[k])` ist der Kern des Algorithmus.

### Komplexität

| Methode | Zeitkomplexität |
|---------|----------------|
| DFT | O(N²) |
| FFT | O(N log N) |

Für N = 1024: FFT ist ~100× schneller als DFT.

### Voraussetzung
N muss eine **Zweierpotenz** sein. Für andere N: automatischer Fallback auf DFT oder Zero-Padding via `fft_padded()`.

### Historische Anmerkung
Cooley und Tukey entdeckten den Algorithmus 1965 neu. Gauß kannte eine äquivalente Methode bereits ca. 1805.

### `ifft(x) → list`

**Inverse FFT** via Konjugationstrick (gleich wie IDFT, aber mit FFT statt DFT):
```
IFFT(X) = konj(FFT(konj(X))) / N
```

### `fft_padded(x) → list`

FFT mit **automatischem Zero-Padding** auf die nächste Zweierpotenz ≥ N.

Zero-Padding verbessert die **Frequenzauflösung** im Spektrum.

---

## Fourier-Reihen

### `fourier_coefficients(f, period, n_terms, n_samples=1000) → (a_coeffs, b_coeffs)`

Berechnet Fourier-Koeffizienten einer T-periodischen Funktion:

**Fourier-Reihe:**
```
f(t) = a₀/2 + Σ_{k=1}^N [aₖ·cos(2πkt/T) + bₖ·sin(2πkt/T)]
```

**Koeffizienten (numerisch via Trapezregel):**
```
aₖ = (2/T) ∫_0^T f(t)·cos(2πkt/T) dt
bₖ = (2/T) ∫_0^T f(t)·sin(2πkt/T) dt
```

### `fourier_series_eval(t, a_coeffs, b_coeffs, period) → float`

Wertet die Fourier-Reihe an einem Punkt `t` aus:
```
f(t) ≈ a₀/2 + Σ_{k=1}^N [aₖ·cos(2πkt/T) + bₖ·sin(2πkt/T)]
```

---

## Fensterfunktionen

Fensterfunktionen reduzieren **Leck-Effekte** (spectral leakage), die entstehen wenn das Signal nicht periodisch mit der FFT-Länge ist.

### `window_hanning(n) → list`

```
w[k] = 0.5·(1 - cos(2πk/N))
```

**Eigenschaften:** Guter Kompromiss, weit verbreitet. Seitenkeulendämpfung: 31 dB.

### `window_hamming(n) → list`

```
w[k] = 0.54 - 0.46·cos(2πk/N)
```

**Eigenschaften:** Leicht bessere Seitenkeulendämpfung als Hanning: 42 dB.

### `window_blackman(n) → list`

```
w[k] = 0.42 - 0.5·cos(2πk/N) + 0.08·cos(4πk/N)
```

**Eigenschaften:** Sehr gute Seitenkeulendämpfung: 74 dB, aber breitere Hauptkeule.

### Fenstervergleich

| Fenster | Hauptkeulenbreite | Seitenkeulendämpfung |
|---------|------------------|---------------------|
| Rechteck | Schmal | 13 dB (schlecht) |
| Hanning | Mittel | 31 dB |
| Hamming | Mittel | 42 dB |
| Blackman | Breit | 74 dB |

### `apply_window(x, window) → list`

Multipliziert das Signal komponentenweise mit der Fensterfunktion:
```
x_windowed[k] = x[k] · window[k]
```

---

## Spektralanalyse

### `power_spectrum(x, sample_rate=1.0) → (freqs, power)`

**Leistungsdichtespektrum (PSD):**
```
P[k] = |X[k]|² / N²
```

Nur das einseitige Spektrum (positive Frequenzen) wird zurückgegeben.

**Frequenzachse:**
```
f[k] = k · sample_rate / N_FFT
```

### `dominant_frequency(x, sample_rate=1.0) → float`

Findet die Frequenz mit der höchsten Leistung (außer DC-Komponente bei k=0).

---

## Kurzzeit-Fourier-Transformation (STFT)

### `stft(x, window_size, hop_size, window_type='hanning') → list`

Analysiert, wie sich das Frequenzspektrum **im Zeitverlauf** ändert:

```
1. Signal in überlappende Fenster aufteilen
2. Jedes Fenster mit Fensterfunktion multiplizieren
3. FFT auf jedes gefensterte Segment anwenden
```

**Parameter:**
| Parameter | Bedeutung |
|-----------|-----------|
| `window_size` | Fensterlänge in Samples |
| `hop_size` | Schrittweite zwischen Fenstern (< window_size → Überlappung) |
| `window_type` | 'hanning', 'hamming', 'blackman', 'rect' |

**Rückgabe:** Liste von FFT-Koeffizienten (Spektrogramm-Zeilen).

---

## Abhängigkeiten

| Modul | Zweck |
|-------|-------|
| `math` | Grundlegende math. Funktionen |
| `cmath` | Komplexe Exponentialfunktionen |
| `numpy` | Nur Hilfsfunktionen (könnte entfernt werden) |
| `typing` | Typ-Annotationen |

---

## Verwendungsbeispiele

```python
from fourier import fft, ifft, fourier_coefficients, power_spectrum, stft
import math

# FFT eines Sinussignals
N = 64
omega = 2 * math.pi * 5 / N   # 5 Hz bei Abtastrate 1
signal = [math.sin(omega * n) for n in range(N)]

X = fft(signal)
# X[5] und X[N-5] sollten die größten Amplituden haben

# Rücktransformation
signal_back = ifft(X)
print([abs(x.real - s) for x, s in zip(signal_back, signal)])  # Nahe 0

# Fourier-Koeffizienten einer Rechteckwelle
def square_wave(t):
    return 1.0 if t % 1.0 < 0.5 else -1.0

a, b = fourier_coefficients(square_wave, period=1.0, n_terms=10)
# b[1], b[3], b[5], ... sollten groß sein (nur ungerade Harmonische)

# Leistungsspektrum
freqs, power = power_spectrum(signal, sample_rate=1.0)

# STFT für zeitvariables Signal
frames = stft(signal, window_size=16, hop_size=8, window_type='hanning')
print(f"Anzahl Frames: {len(frames)}")
```
