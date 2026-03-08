# specialist-maths - Mathematik-Spezialist

## Überblick
Dieses Projekt dient der systematischen Entwicklung mathematischer Expertise.
Es enthält Python-Bibliotheken für alle Gebiete der Mathematik, vollständig
dokumentiert, getestet und als Lernmaterial geeignet.

## Anforderungen
- Python 3.13+
- sympy >= 1.13
- numpy >= 2.2
- scipy >= 1.15
- pytest (für Tests)

## Installation
```bash
pip install sympy numpy scipy matplotlib pytest --break-system-packages
```

## Verfügbare Module

### `src/algebra.py` - Algebra
```python
from algebra import Polynomial, solve_quadratic, gcd, is_prime, euler_phi

# Polynom x^2 - 5x + 6
p = Polynomial([1, -5, 6])
print(p.evaluate(2))   # 0 (Nullstelle)
print(p.derivative())  # 2x - 5

# Quadratische Gleichung lösen
roots = solve_quadratic(1, -5, 6)  # [3.0, 2.0]

# Zahlentheorie
print(gcd(48, 18))           # 6
print(euler_phi(12))          # 4
print(is_prime(97))           # True
```

### `src/analysis.py` - Analysis
```python
import math
from analysis import numerical_derivative, numerical_integral, newton_raphson, taylor_series

# Ableitung von sin(x) bei x=0
print(numerical_derivative(math.sin, 0))   # ≈ 1.0 (= cos(0))

# Integral ∫₀^π sin(x) dx = 2
print(numerical_integral(math.sin, 0, math.pi))  # ≈ 2.0

# Nullstelle von x² - 4
root = newton_raphson(lambda x: x**2 - 4, x0=3)  # ≈ 2.0

# Taylor-Reihe von e^x bei x=1 (Grad 15)
approx = taylor_series(math.exp, center=0, degree=15, evaluate_at=1)  # ≈ e
```

### `src/linear_algebra.py` - Lineare Algebra
```python
from linear_algebra import Matrix, Vector, gram_schmidt

# Vektoren
v1 = Vector([1, 2, 3])
v2 = Vector([4, 5, 6])
print(v1.dot(v2))   # 32
print(v1.norm())    # sqrt(14)

# Matrizen
A = Matrix([[1, 2], [3, 4]])
print(A.determinant())  # -2
b = Vector([5, 6])
x = A.solve(b)          # Löst Ax = b

# Eigenwerte
eigenvalues = A.eigenvalues()

# Gram-Schmidt-Orthogonalisierung
vecs = [Vector([1,1,0]), Vector([1,0,1])]
orthonormal = gram_schmidt(vecs, normalize=True)
```

## Tests ausführen
```bash
python3 -m pytest tests/ -v
```

## Projektstruktur
```
specialist-maths/
├── src/                          # Quellcode
│   ├── algebra.py                # Polynome, Gleichungen, Zahlentheorie
│   ├── analysis.py               # Analysis, Differentiation, Integration
│   ├── linear_algebra.py         # Lineare Algebra, LU/QR/SVD
│   ├── statistics_math.py        # Statistik, Verteilungen, Tests
│   ├── ode.py                    # Differentialgleichungen, Laplace
│   ├── complex_analysis.py       # Komplexe Analysis, Riemann-Zeta
│   ├── analytic_number_theory.py # Analytische Zahlentheorie
│   ├── proof_theory.py           # Beweistheorie, Millennium-Probleme
│   ├── fourier.py                # Fourier-Transformation (FFT)
│   ├── numerical_methods.py      # Interpolation, Optimierung, Simplex
│   └── build.txt                 # Build-Nummer (aktuell: 8)
├── tests/                        # 333 Tests (100% grün)
├── research/                     # Mathematische Notizen und Recherche
├── dev-log/                      # Entwicklungs-Tagebuch
├── debugging/                    # Debug-Skripte
├── build/                        # Ausführbare Dateien
├── FEATURES.md                   # Feature-Liste
├── BUGS.md                       # Bug-Tracker
└── OPTIMIZE.md                   # Optimierungsideen
```

## Mathematische Themengebiete
1. **Algebra**: Polynomrechnung, lineare/quadratische Gleichungen, Zahlentheorie
2. **Analysis**: Differentiation, Integration, Taylor-Reihen, Nullstellensuche
3. **Lineare Algebra**: Matrizen, Vektoren, Eigenwerte, LU/QR/SVD
4. **Statistik**: Verteilungen, Hypothesentests, Monte-Carlo
5. **Differentialgleichungen**: Euler, RK4/RK45, Laplace-Transformation
6. **Komplexe Analysis**: Riemann-Zeta, Gamma-Funktion, Nullstellen
7. **Analytische Zahlentheorie**: π(x), Li(x), Primzahlsatz, Dirichlet-L-Reihen
8. **Beweistheorie**: Collatz, Goldbach, Millennium-Probleme
9. **Fourier-Analysis**: DFT/FFT, Fourier-Reihen, STFT, Fensterfunktionen
10. **Numerische Methoden**: Interpolation (Lagrange/Newton/Splines), BFGS, Simplex
11. **Kryptographie**: RSA (Schlüsselerzeugung, Ver-/Entschlüsselung)

---
*Autor: Kurt Ingwer*
*Build: 8 | Stand: 2026-03-08 | Tests: 333/333 grün*
