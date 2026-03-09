# latex_export.py — Dokumentation

## Übersicht

Das Modul `latex_export.py` stellt Funktionen bereit, um mathematische Ergebnisse
in LaTeX-Code zu exportieren. Der erzeugte LaTeX-Code ist direkt kompilierbar
mit `pdflatex`, `lualatex` oder `xelatex`.

## Abhängigkeiten

- `sympy` — für `sympy_to_latex()`
- `math` — für Spezialwert-Erkennung (π, e, √n)
- Python-Standardbibliothek (`math.gcd`)

---

## Funktionsreferenz

### `number_to_latex(x, precision=6) -> str`

Konvertiert eine Zahl in einen LaTeX-String.

| Eingabe | Ausgabe |
|---------|---------|
| `42` | `"42"` |
| `-7` | `"-7"` |
| `math.pi` | `"\pi"` |
| `math.sqrt(2)` | `"\sqrt{2}"` |
| `0.5` | `"\frac{1}{2}"` |
| `3+4j` | `"3 + 4i"` |

**Erkannte Spezialwerte:** $\pi$, $2\pi$, $\frac{\pi}{2}$, $\frac{\pi}{3}$, $\frac{\pi}{4}$, $e$, $\sqrt{2}$, $\sqrt{3}$, $\sqrt{5}$, $\sqrt{6}$, $\sqrt{7}$

---

### `fraction_to_latex(num, den) -> str`

Erstellt einen LaTeX-Bruch $\frac{\text{num}}{\text{den}}$.

Vereinfacht automatisch durch den ggT. Normalisiert das Vorzeichen in den Zähler.

```python
fraction_to_latex(1, 2)   # "\frac{1}{2}"
fraction_to_latex(2, 4)   # "\frac{1}{2}"  (vereinfacht)
fraction_to_latex(6, 3)   # "2"            (Ganzzahl)
fraction_to_latex(1, -3)  # "\frac{-1}{3}" (Vorzeichen normalisiert)
```

**Formel:** $\frac{p}{q}$ mit $\gcd(|p|, q) = 1$, $q > 0$

---

### `polynomial_to_latex(coefficients, var='x') -> str`

Konvertiert eine Koeffizientenliste in LaTeX-Polynom-Darstellung.

**Konvention:** `coefficients[0]` ist der höchste Grad.

```python
polynomial_to_latex([1, -3, 2])     # "x^{2} - 3x + 2"
polynomial_to_latex([2, 0, -1])     # "2x^{2} - 1"
polynomial_to_latex([1, 0, 0], 't') # "t^{2}"
```

**Ausgabe:** $a_n x^n + a_{n-1} x^{n-1} + \ldots + a_1 x + a_0$

Sonderregeln:
- Koeffizient $\pm 1$ bei Termen mit Variable: wird weggelassen (`x^{2}` statt `1x^{2}`)
- Koeffizient $0$: Term wird weggelassen
- Exponent $1$: `x` statt `x^{1}`
- Exponent $0$: nur Konstante

---

### `matrix_to_latex(matrix_data, env='pmatrix') -> str`

Konvertiert eine 2D-Liste in eine LaTeX-Matrix.

| Umgebung | Klammern |
|----------|----------|
| `pmatrix` | Rund ( ) |
| `bmatrix` | Eckig [ ] |
| `vmatrix` | Betragsstriche \| \| |
| `Bmatrix` | Geschwungen { } |
| `matrix` | Keine |

```python
matrix_to_latex([[1, 2], [3, 4]])
# \begin{pmatrix}
# 1 & 2
# 3 & 4
# \end{pmatrix}
```

**LaTeX-Ausgabe:**
$$\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$$

---

### `vector_to_latex(components, column=True) -> str`

```python
vector_to_latex([1, 2, 3])             # Spaltenvektor
vector_to_latex([1, 2, 3], column=False)  # Zeilenvektor
```

**Spaltenvektor:** $\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$

**Zeilenvektor:** $\begin{pmatrix} 1 & 2 & 3 \end{pmatrix}$

---

### `equation_to_latex(lhs, rhs, eq_type='=') -> str`

Erstellt LaTeX-Gleichungen verschiedener Typen.

| `eq_type` | Symbol | LaTeX |
|-----------|--------|-------|
| `'='` | $=$ | `=` |
| `'<='` | $\leq$ | `\leq` |
| `'>='` | $\geq$ | `\geq` |
| `'!='` | $\neq$ | `\neq` |
| `'~='` | $\approx$ | `\approx` |
| `'iff'` | $\iff$ | `\iff` |
| `'implies'` | $\implies$ | `\implies` |

---

### `sympy_to_latex(expr) -> str`

Delegiert an `sympy.latex()` für vollständige SymPy-Ausdrücke.

```python
import sympy as sp
x = sp.Symbol('x')
sympy_to_latex(sp.integrate(sp.sin(x), x))  # "-\cos{x}"
```

---

### `integral_to_latex(integrand, var, lower=None, upper=None) -> str`

```python
integral_to_latex("x^2", "x", 0, 1)   # "\int_{0}^{1} x^2 \, dx"
integral_to_latex("x^2", "x")          # "\int x^2 \, dx"
```

**Bestimmtes Integral:** $\int_{a}^{b} f(x) \, dx$

**Unbestimmtes Integral:** $\int f(x) \, dx$

---

### `sum_to_latex(term, index, lower, upper) -> str`

```python
sum_to_latex("n^2", "n", "1", r"\infty")
# "\sum_{n=1}^{\infty} n^2"
```

**Ausgabe:** $\sum_{n=1}^{\infty} n^2$

---

### `limit_to_latex(expr, var, point, direction='') -> str`

```python
limit_to_latex("sin(x)/x", "x", "0")     # \lim_{x \to 0} sin(x)/x
limit_to_latex("f(x)", "x", "0", "+")    # \lim_{x \to 0^{+}} f(x)
limit_to_latex("f(x)", "x", "0", "-")    # \lim_{x \to 0^{-}} f(x)
```

**Zweiseitig:** $\lim_{x \to 0} \frac{\sin(x)}{x}$

**Einseitig:** $\lim_{x \to 0^+} f(x)$

---

### `derivative_to_latex(func, var, order=1) -> str`

```python
derivative_to_latex("f(x)", "x")     # "\frac{d}{dx} f(x)"
derivative_to_latex("f(x)", "x", 3)  # "\frac{d^{3}}{dx^{3}} f(x)"
```

**Erste Ableitung:** $\frac{d}{dx} f(x)$

**n-te Ableitung:** $\frac{d^n}{dx^n} f(x)$

---

### `theorem_to_latex(name, statement, proof='') -> str`

```python
theorem_to_latex("Pythagoras",
    r"a^2 + b^2 = c^2",
    r"Im rechtwinkligen Dreieck...")
```

**Ausgabe:**
```latex
\begin{theorem}[Pythagoras]
a^2 + b^2 = c^2
\end{theorem}
\begin{proof}
Im rechtwinkligen Dreieck...
\end{proof}
```

> Benötigt `\usepackage{amsthm}` in der Präambel.

---

### `solution_to_latex(result, title='') -> str`

Konvertiert ein Dictionary in eine `align*`-Umgebung:

```python
solution_to_latex({"x": 3, "y": -1}, "Lösung")
```

**Ausgabe:**
```latex
\textbf{Lösung}

\begin{align*}
  x &= 3 \\
  y &= -1
\end{align*}
```

---

### `full_document(content, title, author, packages) -> str`

Erstellt ein vollständiges LaTeX-Dokument:

```python
doc = full_document(
    r"\section{Algebra}\n...",
    title="Mathematische Notizen",
    author="Kurt Ingwer",
    packages=["graphicx"]
)
save_latex(doc, "ausgabe.tex")
```

**Standard-Pakete:** `amsmath`, `amsfonts`, `amssymb`, `amsthm`, `geometry`

---

### `save_latex(content, filepath) -> None`

Speichert einen LaTeX-String in eine Datei (UTF-8).

```python
save_latex(full_document(content), "/tmp/math.tex")
```

---

## Vollständiges Beispiel

```python
import latex_export as lx
import math

# Polynom
poly = lx.polynomial_to_latex([1, -3, 2])          # x^{2} - 3x + 2

# Integral
intg = lx.integral_to_latex("x^2 - 3x + 2", "x", 0, 2)

# Theorem
thm = lx.theorem_to_latex(
    "Nullstellensatz",
    f"Das Polynom ${poly}$ hat Nullstellen bei $x_1=1, x_2=2$."
)

# Vollständiges Dokument
doc = lx.full_document(
    f"\\section{{Analyse}}\nDas Integral:\n$${intg}$$\n\n{thm}",
    title="Algebraische Analysis",
    author="Kurt Ingwer"
)

lx.save_latex(doc, "ergebnis.tex")
```

---

*Letzte Änderung: 2026-03-09 | Autor: Kurt Ingwer*
