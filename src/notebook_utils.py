"""
@file notebook_utils.py
@brief Jupyter-Notebook-kompatible Ausgabe-Hilfsfunktionen.
@description
    Dieses Modul stellt Hilfsfunktionen bereit, um mathematische Ergebnisse
    in ansprechender Form in Jupyter-Notebooks und interaktiven Umgebungen
    darzustellen. Es unterstützt HTML-Formatierung, LaTeX-Ausgabe und die
    Erzeugung von Demo-Notebooks.

    ## Verwendung in Jupyter

    ```python
    from notebook_utils import display_math, matrix_to_html
    import sympy as sp

    x = sp.Symbol('x')
    display_math(sp.integrate(sp.sin(x), x))  # Zeigt -cos(x)

    import numpy as np
    A = np.array([[1, 2], [3, 4]])
    print(matrix_to_html(A))  # HTML-Tabelle
    ```

    ## Fallback im Terminal

    Alle Funktionen funktionieren auch ohne Jupyter – sie geben dann
    eine textuelle Darstellung aus.

@author Michael Fuhrmann
@lastModified 2026-03-10
@version 1.0.0
"""

import json
import os
from typing import Any, Dict, List, Optional, Union


# ===========================================================================
# HTML-Formatierungsfunktionen
# ===========================================================================

def matrix_to_html(matrix) -> str:
    """
    @brief Formatiert eine Matrix als HTML-Tabelle (Jupyter-ready).
    @description
        Erzeugt eine HTML-Tabelle, die eine Matrix (numpy-Array, Liste von
        Listen oder SymPy-Matrix) darstellt. Der Stil ist für Jupyter-Notebooks
        optimiert.

        Das Ergebnis kann direkt mit ``IPython.display.HTML()`` angezeigt werden.

    @param matrix  Zu formatierende Matrix (numpy-Array, Liste oder SymPy-Matrix).
    @return        HTML-String der Matrix-Tabelle.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Matrix in eine Liste von Listen umwandeln
    rows = _matrix_to_rows(matrix)

    if not rows:
        return "<table><tr><td>(leere Matrix)</td></tr></table>"

    # HTML-Tabelle mit Jupyter-Stil aufbauen
    html_parts: List[str] = []
    html_parts.append(
        '<table style="border-collapse: collapse; font-family: monospace;">'
    )

    for row in rows:
        html_parts.append('  <tr>')
        for cell in row:
            # Zellwert formatieren (Fließkommazahlen auf 4 Stellen kürzen)
            formatted = _format_cell(cell)
            html_parts.append(
                f'    <td style="border: 1px solid #ccc; padding: 4px 8px; '
                f'text-align: right;">{formatted}</td>'
            )
        html_parts.append('  </tr>')

    html_parts.append('</table>')
    return '\n'.join(html_parts)


def _matrix_to_rows(matrix) -> List[List[Any]]:
    """
    @brief Konvertiert verschiedene Matrix-Formate in eine Liste von Zeilen.
    @param matrix  Eingabematrix (numpy, Liste oder SymPy).
    @return        Liste von Zeilen, jede Zeile eine Liste von Werten.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    try:
        import numpy as np
        if isinstance(matrix, np.ndarray):
            if matrix.ndim == 1:
                # 1D-Array als einzeilige Matrix darstellen
                return [[matrix.tolist()]]
            elif matrix.ndim == 2:
                return matrix.tolist()
            else:
                # Höherdimensional: erste 2D-Schicht nehmen
                return matrix.reshape(-1, matrix.shape[-1]).tolist()
    except ImportError:
        pass

    # Liste von Listen
    if isinstance(matrix, list):
        if not matrix:
            return []
        if isinstance(matrix[0], list):
            return matrix
        else:
            return [matrix]  # Flache Liste → einzeilig

    # SymPy-Matrix
    try:
        import sympy as sp
        if isinstance(matrix, sp.Matrix):
            rows_s = matrix.rows
            cols_s = matrix.cols
            return [[matrix[i, j] for j in range(cols_s)]
                    for i in range(rows_s)]
    except ImportError:
        pass

    # Fallback
    return [[str(matrix)]]


def _format_cell(value: Any) -> str:
    """
    @brief Formatiert einen einzelnen Zellwert für HTML.
    @param value  Zellwert (Zahl, SymPy-Ausdruck, etc.).
    @return       Formatierter String.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    if isinstance(value, float):
        # Fließkommazahl: 6 signifikante Stellen
        return f"{value:.6g}"
    elif isinstance(value, complex):
        # Komplexe Zahl
        if value.imag == 0:
            return f"{value.real:.6g}"
        return f"{value.real:.4g}+{value.imag:.4g}j"
    else:
        return str(value)


def polynomial_to_latex(coeffs: List[float], variable: str = 'x') -> str:
    """
    @brief Gibt ein Polynom als LaTeX-String zurück.
    @description
        Erzeugt aus einer Koeffizientenliste (aufsteigend nach Grad)
        einen LaTeX-Ausdruck für das Polynom.

        Konvention: coeffs[k] ist der Koeffizient von x^k.
        Beispiel: coeffs=[1, 0, -2, 3] → 3x³ - 2x² + 1

    @param coeffs    Liste der Koeffizienten [a_0, a_1, ..., a_n].
    @param variable  Name der Variablen (Standard: 'x').
    @return          LaTeX-String des Polynoms.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    if not coeffs:
        return "0"

    # Terme von höchstem Grad nach unten aufbauen
    terms: List[str] = []
    degree = len(coeffs) - 1

    for k in range(degree, -1, -1):
        coeff = coeffs[k]

        # Null-Terme überspringen
        if coeff == 0:
            continue

        # Vorzeichen bestimmen
        is_first = (len(terms) == 0)
        if coeff < 0:
            sign = "-" if is_first else " - "
            abs_coeff = abs(coeff)
        else:
            sign = "" if is_first else " + "
            abs_coeff = coeff

        # Koeffizient formatieren
        if isinstance(abs_coeff, float):
            coeff_str = f"{abs_coeff:.6g}" if abs_coeff != 1.0 or k == 0 else ""
        else:
            coeff_str = str(abs_coeff) if abs_coeff != 1 or k == 0 else ""

        # Potenzterm
        if k == 0:
            # Konstantterm
            terms.append(f"{sign}{abs_coeff:.6g}" if isinstance(abs_coeff, float)
                        else f"{sign}{abs_coeff}")
        elif k == 1:
            # Linearer Term
            if coeff_str:
                terms.append(f"{sign}{coeff_str}{variable}")
            else:
                terms.append(f"{sign}{variable}")
        else:
            # Höhere Potenzen
            if coeff_str:
                terms.append(f"{sign}{coeff_str}{variable}^{{{k}}}")
            else:
                terms.append(f"{sign}{variable}^{{{k}}}")

    if not terms:
        return "0"

    return "".join(terms)


def result_to_html(result: Dict[str, Any]) -> str:
    """
    @brief Formatiert ein Berechnungsergebnis als schönes HTML.
    @description
        Erwartet ein Dictionary mit beliebigen Schlüssel-Wert-Paaren und
        erzeugt eine formatierte HTML-Darstellung. Numerische Werte werden
        auf 8 signifikante Stellen gerundet, Listen werden als Aufzählungen
        dargestellt.

    @param result  Dictionary mit Berechnungsergebnissen.
    @return        HTML-String der formatierten Ausgabe.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    html_parts: List[str] = []
    html_parts.append(
        '<div style="font-family: monospace; background: #f8f8f8; '
        'padding: 12px; border-left: 4px solid #2196F3; margin: 8px 0;">'
    )
    html_parts.append('<table style="border-collapse: collapse; width: 100%;">')

    for key, value in result.items():
        # Schlüssel formatieren
        key_html = f'<strong style="color: #1a237e;">{key}</strong>'

        # Wert formatieren
        if isinstance(value, float):
            val_html = f"{value:.8g}"
        elif isinstance(value, (list, tuple)):
            items = ", ".join(_format_cell(v) for v in value)
            val_html = f"[{items}]"
        elif isinstance(value, dict):
            items = ", ".join(f"{k}: {_format_cell(v)}" for k, v in value.items())
            val_html = f"{{{items}}}"
        elif isinstance(value, bool):
            color = "#2e7d32" if value else "#c62828"
            val_html = f'<span style="color: {color};">{value}</span>'
        else:
            val_html = str(value)

        html_parts.append(
            f'  <tr>'
            f'    <td style="padding: 4px 12px 4px 4px; vertical-align: top; '
            f'color: #555; width: 40%;">{key_html}</td>'
            f'    <td style="padding: 4px;">{val_html}</td>'
            f'  </tr>'
        )

    html_parts.append('</table>')
    html_parts.append('</div>')
    return '\n'.join(html_parts)


def display_math(expr) -> None:
    """
    @brief Zeigt einen mathematischen Ausdruck in Notebook oder Terminal an.
    @description
        Versucht zunächst, ``IPython.display.Math`` zu verwenden, um den
        Ausdruck als LaTeX-Formel in einem Jupyter-Notebook anzuzeigen.
        Falls IPython nicht verfügbar ist (Terminal-Betrieb), wird der
        Ausdruck als String ausgegeben.

        Für SymPy-Ausdrücke wird ``sympy.latex()`` zur LaTeX-Konvertierung
        verwendet.

    @param expr  Mathematischer Ausdruck (SymPy, Zahl, String).
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Versuche LaTeX-Darstellung zu erzeugen
    latex_str = _to_latex_str(expr)

    # Versuche IPython-Anzeige (Jupyter-Notebook)
    try:
        from IPython.display import Math, display  # type: ignore[import]
        display(Math(latex_str))
    except ImportError:
        # Fallback: Terminalausgabe
        print(f"  {latex_str}")
    except Exception:
        # Allgemeiner Fallback
        print(str(expr))


def _to_latex_str(expr) -> str:
    """
    @brief Konvertiert einen Ausdruck in einen LaTeX-String.
    @param expr  Mathematischer Ausdruck.
    @return      LaTeX-String.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # SymPy-Ausdruck
    try:
        import sympy as sp
        if isinstance(expr, (sp.Basic, sp.Matrix)):
            return sp.latex(expr)
    except ImportError:
        pass

    # Zahlen direkt formatieren
    if isinstance(expr, float):
        return f"{expr:.6g}"
    elif isinstance(expr, complex):
        return f"{expr.real:.4g} + {expr.imag:.4g}i"

    # String: bereits LaTeX
    if isinstance(expr, str):
        return expr

    return str(expr)


# ===========================================================================
# Notebook-Erzeugung
# ===========================================================================

def save_notebook_demo(output_path: str = "notebooks/demo.ipynb") -> None:
    """
    @brief Erzeugt ein Demo-Jupyter-Notebook mit Beispielberechnungen.
    @description
        Erstellt ein vollständiges Jupyter-Notebook (.ipynb) im nbformat 4.4,
        das die wichtigsten specialist-maths Module demonstriert:
        - Algebra: Polynome und Primzahlen
        - Analysis: Differentiation und Integration
        - Lineare Algebra: Eigenwerte
        - Visualisierung: Einfache Plots

        Das Notebook kann direkt mit ``jupyter notebook`` geöffnet werden.

    @param output_path  Pfad zur Ausgabedatei (Standard: "notebooks/demo.ipynb").
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    # Ausgabeverzeichnis erstellen
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Notebook-Struktur gemäß nbformat 4.4
    notebook = {
        "nbformat": 4,
        "nbformat_minor": 4,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "name": "python",
                "version": "3.13.0"
            }
        },
        "cells": [
            _make_markdown_cell("# specialist-maths – Demo-Notebook\n\n"
                "Dieses Notebook demonstriert die wichtigsten Funktionen des "
                "specialist-maths Projekts.\n\n"
                "**Autor:** Michael Fuhrmann  \n"
                "**Stand:** 2026-03-10"),

            _make_markdown_cell("## 1. Setup: Pfade und Imports"),

            _make_code_cell(
                "import sys\n"
                "sys.path.insert(0, '../src')\n\n"
                "# Alle Kernmodule importieren\n"
                "from algebra import Polynomial, is_prime, prime_factorization\n"
                "from analysis import numerical_derivative, numerical_integral\n"
                "from linear_algebra import Matrix, Vector\n"
                "from notebook_utils import matrix_to_html, polynomial_to_latex, "
                "result_to_html, display_math\n\n"
                "print('Imports erfolgreich!')"
            ),

            _make_markdown_cell("## 2. Algebra: Polynome und Primzahlen"),

            _make_code_cell(
                "# Polynom p(x) = x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)\n"
                "p = Polynomial([1, -6, 11, -6])  # Koeffizienten aufsteigend\n"
                "latex_str = polynomial_to_latex([-6, 11, -6, 1])\n"
                "print(f'LaTeX: $p(x) = {latex_str}$')\n\n"
                "# Primzahlen überprüfen\n"
                "primes = [n for n in range(2, 50) if is_prime(n)]\n"
                "print(f'Primzahlen bis 50: {primes}')\n\n"
                "# Primfaktorzerlegung\n"
                "print(f'Primfaktoren von 360: {prime_factorization(360)}')"
            ),

            _make_markdown_cell("## 3. Analysis: Differentiation und Integration"),

            _make_code_cell(
                "import math\n\n"
                "# Numerische Ableitung von sin(x) bei x=0 (Erwartung: 1.0)\n"
                "df = numerical_derivative(math.sin, 0.0)\n"
                "print(f\"sin'(0) ≈ {df:.8f}  (exakt: 1.0)\")\n\n"
                "# Numerische Integration von x^2 von 0 bis 1 (Erwartung: 1/3)\n"
                "integral = numerical_integral(lambda x: x**2, 0.0, 1.0)\n"
                "print(f'∫x² dx [0,1] ≈ {integral:.8f}  (exakt: 0.33333...)')\n\n"
                "# Ergebnis als HTML-Tabelle\n"
                "result = {\n"
                "    \"sin'(0)\": df,\n"
                "    \"∫x² [0,1]\": integral,\n"
                "    \"Fehler sin'\": abs(df - 1.0),\n"
                "    \"Fehler ∫\": abs(integral - 1/3)\n"
                "}\n"
                "print(result_to_html(result))"
            ),

            _make_markdown_cell("## 4. Lineare Algebra: Matrizen und Eigenwerte"),

            _make_code_cell(
                "import numpy as np\n\n"
                "# 3x3-Hilbert-Matrix\n"
                "H = Matrix([[1, 1/2, 1/3],\n"
                "            [1/2, 1/3, 1/4],\n"
                "            [1/3, 1/4, 1/5]])\n\n"
                "print('Hilbert-Matrix H_3:')\n"
                "print(matrix_to_html(H.data))\n\n"
                "print(f'Determinante: {H.determinant():.6e}')\n"
                "print(f'Eigenwerte: {[f\"{e:.4f}\" for e in H.eigenvalues()]}')"
            ),

            _make_markdown_cell("## 5. SymPy: Symbolische Mathematik"),

            _make_code_cell(
                "import sympy as sp\n"
                "from notebook_utils import display_math\n\n"
                "x = sp.Symbol('x')\n\n"
                "# Integral berechnen\n"
                "integral = sp.integrate(sp.sin(x)**2, x)\n"
                "print('∫sin²(x) dx =')\n"
                "display_math(integral)\n\n"
                "# Grenzwert\n"
                "lim = sp.limit(sp.sin(x)/x, x, 0)\n"
                "print(f'lim(x→0) sin(x)/x = {lim}')\n\n"
                "# Taylorentwicklung\n"
                "taylor = sp.series(sp.exp(x), x, 0, 6)\n"
                "print('Taylor e^x um 0 (6 Terme):')\n"
                "display_math(taylor)"
            ),

            _make_markdown_cell("## 6. Plugin-System"),

            _make_code_cell(
                "from plugin_system import discover_plugins, get_plugin, list_plugins\n\n"
                "# Plugins laden\n"
                "loaded = discover_plugins('../plugins/')\n"
                "print(f'Geladene Plugins: {loaded}')\n\n"
                "# Beispiel-Plugin nutzen\n"
                "if 'example_plugin' in list_plugins():\n"
                "    plugin = get_plugin('example_plugin')\n"
                "    print(f'28 ist vollkommen: {plugin.perfect_number_check(28)}')\n"
                "    print(f'Quersumme(1234) = {plugin.digit_sum(1234)}')\n"
                "    print(f'Collatz-Länge(27) = {plugin.collatz_length(27)}')"
            ),

            _make_markdown_cell(
                "---\n*Notebook automatisch generiert von specialist-maths v31*"
            ),
        ]
    }

    # Notebook als JSON-Datei speichern
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(notebook, f, indent=2, ensure_ascii=False)

    print(f"Demo-Notebook erstellt: {output_path}")


def _make_code_cell(source: str) -> Dict[str, Any]:
    """
    @brief Erzeugt eine Code-Zelle für ein Jupyter-Notebook.
    @param source  Python-Quellcode der Zelle.
    @return        nbformat-konformes Dictionary.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source
    }


def display_matrix_html(M) -> str:
    """
    @brief Zeigt eine Matrix als HTML-Tabelle an (Jupyter-kompatibel).
    @description
        Konvertiert eine Matrix (numpy-Array, Liste oder SymPy-Matrix) in
        eine HTML-Tabelle und zeigt sie via IPython.display an.
        Falls IPython nicht verfügbar ist, wird der HTML-String zurückgegeben.

        Unterschied zu matrix_to_html(): Diese Funktion gibt das Ergebnis
        direkt an display() weiter und gibt bei Jupyter nichts zurück,
        bei Terminal-Betrieb den HTML-String.

    @param M   Matrix (numpy-Array, Liste von Listen oder SymPy-Matrix)
    @return    HTML-String (nur im Terminal-Modus), None in Jupyter
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # HTML-Tabelle generieren
    html_str = matrix_to_html(M)

    # Versuche IPython-Anzeige (Jupyter-Notebook)
    try:
        from IPython.display import HTML, display  # type: ignore[import]
        display(HTML(html_str))
        return html_str  # Im Jupyter-Modus trotzdem zurückgeben (für Tests)
    except ImportError:
        # Terminal-Modus: HTML-String zurückgeben
        return html_str
    except Exception:
        # Allgemeiner Fallback
        return html_str


def display_polynomial_latex(coeffs: List[float], variable: str = 'x') -> str:
    """
    @brief Zeigt ein Polynom als LaTeX-Formel an (Jupyter-kompatibel).
    @description
        Konvertiert eine Koeffizientenliste in einen LaTeX-String und
        zeigt ihn via IPython.display.Math an.
        Falls IPython nicht verfügbar ist, wird der LaTeX-String zurückgegeben.

        Koeffizientenkonvention: coeffs[k] ist der Koeffizient von x^k.
        Beispiel: coeffs=[6, -5, 1] → x² - 5x + 6

    @param coeffs    Liste der Koeffizienten [a_0, a_1, ..., a_n]
    @param variable  Name der Variablen (Standard: 'x')
    @return          LaTeX-String des Polynoms
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # LaTeX-String aus Koeffizienten erstellen
    latex_str = polynomial_to_latex(coeffs, variable)

    # Versuche IPython-Anzeige (Jupyter-Notebook)
    try:
        from IPython.display import Math, display  # type: ignore[import]
        display(Math(latex_str))
        return latex_str  # String trotzdem zurückgeben (für Tests)
    except ImportError:
        # Terminal-Modus: LaTeX-String zurückgeben
        return latex_str
    except Exception:
        # Allgemeiner Fallback
        return latex_str


def _make_markdown_cell(source: str) -> Dict[str, Any]:
    """
    @brief Erzeugt eine Markdown-Zelle für ein Jupyter-Notebook.
    @param source  Markdown-Text der Zelle.
    @return        nbformat-konformes Dictionary.
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": source
    }
