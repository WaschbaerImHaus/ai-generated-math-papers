"""
@file latex_export.py
@brief LaTeX-Export-Modul für mathematische Ergebnisse.
@description
    Ermöglicht den Export mathematischer Berechnungsergebnisse als LaTeX-Code.
    Unterstützt Zahlen, Polynome, Matrizen, Gleichungssysteme, Formeln und
    vollständige LaTeX-Dokumente.

    Alle erzeugten LaTeX-Strings sind direkt in .tex-Dateien verwendbar und
    mit gängigen LaTeX-Compilern (pdflatex, lualatex, xelatex) kompilierbar.

@author Michael Fuhrmann
@date 2026-03-09
@version 1.0.0
"""

import math
import sympy as sp
from typing import Union, Optional, Any
from math import gcd


# =============================================================================
# HILFSFUNKTIONEN
# =============================================================================

def _detect_special_value(x: float) -> Optional[str]:
    """
    @brief Erkennt bekannte mathematische Konstanten und gibt LaTeX zurück.
    @description
        Prüft ob der Wert eine bekannte Konstante ist:
        - π (Pi)
        - e (Eulersche Zahl)
        - √2, √3 usw.
    @param x Zu prüfender Fließkommawert.
    @return LaTeX-String falls erkannt, sonst None.
    @date 2026-03-09
    """
    # Pi und Vielfache
    if abs(x - math.pi) < 1e-10:
        return r"\pi"
    if abs(x + math.pi) < 1e-10:
        return r"-\pi"
    if abs(x - 2 * math.pi) < 1e-10:
        return r"2\pi"
    if abs(x - math.pi / 2) < 1e-10:
        return r"\frac{\pi}{2}"
    if abs(x - math.pi / 3) < 1e-10:
        return r"\frac{\pi}{3}"
    if abs(x - math.pi / 4) < 1e-10:
        return r"\frac{\pi}{4}"

    # Eulersche Zahl
    if abs(x - math.e) < 1e-10:
        return r"e"

    # Häufige Wurzeln
    for n in [2, 3, 5, 6, 7]:
        if abs(x - math.sqrt(n)) < 1e-10:
            return r"\sqrt{" + str(n) + r"}"
        if abs(x + math.sqrt(n)) < 1e-10:
            return r"-\sqrt{" + str(n) + r"}"

    return None


def _float_as_fraction(x: float, max_denom: int = 100) -> Optional[tuple]:
    """
    @brief Versucht einen Float als einfachen Bruch darzustellen.
    @description
        Nutzt den Kettenbruchalgorithmus um einen Fließkommawert als
        p/q mit kleinen Nenner zu erkennen. Gibt None zurück wenn kein
        einfacher Bruch gefunden wurde.
    @param x Zu prüfender Fließkommawert.
    @param max_denom Maximaler Nenner für die Erkennung.
    @return Tupel (Zähler, Nenner) oder None.
    @date 2026-03-09
    """
    # Spezialfall: ganzzahlig
    if abs(x - round(x)) < 1e-10:
        return None

    # Suche einfachen Bruch
    for denom in range(2, max_denom + 1):
        numer = x * denom
        if abs(numer - round(numer)) < 1e-8:
            n = round(numer)
            g = gcd(abs(n), denom)
            return (n // g, denom // g)
    return None


# =============================================================================
# ÖFFENTLICHE FUNKTIONEN
# =============================================================================

def number_to_latex(x: Union[int, float, complex], precision: int = 6) -> str:
    """
    @brief Konvertiert eine Zahl in einen LaTeX-String.
    @description
        Behandelt folgende Fälle:
        - int: Direkte Darstellung ("42", "-7")
        - float: Fließkommawert mit precision Nachkommastellen,
                 oder \frac{} falls einfacher Bruch erkannt,
                 oder Sonderzeichen falls π, e, √n erkannt
        - complex: "a + bi" Darstellung
        - Spezialwerte: π → \pi, e → e, √2 → \sqrt{2}
    @param x Zu konvertierende Zahl.
    @param precision Anzahl Nachkommastellen für Fließkommazahlen.
    @return LaTeX-String der Zahl.
    @date 2026-03-09
    """
    # Komplexe Zahlen
    if isinstance(x, complex):
        real_part = x.real
        imag_part = x.imag

        # Realteil formatieren
        if abs(real_part - round(real_part)) < 1e-10:
            real_str = str(int(round(real_part)))
        else:
            real_str = f"{real_part:.{precision}g}"

        # Imaginärteil formatieren
        if abs(imag_part - round(imag_part)) < 1e-10:
            imag_val = int(round(imag_part))
        else:
            imag_val = imag_part

        if imag_part == 0:
            return real_str
        elif real_part == 0:
            if imag_val == 1:
                return "i"
            elif imag_val == -1:
                return "-i"
            return f"{imag_val}i"
        else:
            if imag_val == 1:
                return f"{real_str} + i"
            elif imag_val == -1:
                return f"{real_str} - i"
            elif isinstance(imag_val, int) and imag_val < 0:
                return f"{real_str} - {-imag_val}i"
            else:
                return f"{real_str} + {imag_val}i"

    # Ganzzahlen
    if isinstance(x, int):
        return str(x)

    # Float: zunächst auf Ganzzahligkeit prüfen
    if isinstance(x, float):
        if abs(x - round(x)) < 1e-10:
            return str(int(round(x)))

        # Spezialwerte prüfen
        special = _detect_special_value(x)
        if special is not None:
            return special

        # Als Bruch darstellen falls möglich
        frac = _float_as_fraction(x)
        if frac is not None:
            return fraction_to_latex(frac[0], frac[1])

        # Standarddarstellung
        return f"{x:.{precision}f}".rstrip('0').rstrip('.')

    # Fallback für andere numerische Typen
    return str(x)


def fraction_to_latex(num: int, den: int) -> str:
    """
    @brief Konvertiert einen Bruch num/den in LaTeX.
    @description
        Vereinfacht den Bruch durch den ggT und gibt \frac{num}{den} zurück.
        Bei den == 1 wird nur der Zähler zurückgegeben.
        Bei negativem Nenner wird das Vorzeichen in den Zähler verschoben.
    @param num Zähler des Bruchs.
    @param den Nenner des Bruchs (darf nicht 0 sein).
    @return LaTeX-String \frac{num}{den} oder vereinfacht.
    @date 2026-03-09
    """
    if den == 0:
        raise ValueError("Nenner darf nicht 0 sein")

    # Vorzeichen normalisieren: Nenner immer positiv
    if den < 0:
        num = -num
        den = -den

    # Vereinfachen durch ggT
    g = gcd(abs(num), den)
    num //= g
    den //= g

    # Ganzzahl falls Nenner = 1
    if den == 1:
        return str(num)

    return r"\frac{" + str(num) + r"}{" + str(den) + r"}"


def polynomial_to_latex(coefficients: list, var: str = 'x') -> str:
    """
    @brief Konvertiert eine Koeffizientenliste in einen LaTeX-Polynomial-String.
    @description
        Koeffizienten[0] ist der höchste Grad. Beispiel:
        [1, -3, 2] → "x^{2} - 3x + 2" (für x^2 - 3x + 2)

        Sonderbehandlungen:
        - Koeffizient 1 oder -1: "x^n" statt "1x^n"
        - Koeffizient 0: Term wird weggelassen
        - Exponent 1: "x" statt "x^{1}"
        - Exponent 0: nur Konstante, kein "x^{0}"
    @param coefficients Liste von Koeffizienten (höchster Grad zuerst).
    @param var Name der Variablen (Standard: 'x').
    @return LaTeX-String des Polynoms.
    @date 2026-03-09
    """
    if not coefficients:
        return "0"

    n = len(coefficients) - 1  # Grad des Polynoms
    terms = []

    for i, coeff in enumerate(coefficients):
        exp = n - i  # Exponent für diesen Term

        # Null-Koeffizienten überspringen
        if coeff == 0:
            continue

        # Betrag und Vorzeichen trennen
        abs_coeff = abs(coeff)
        sign = "-" if coeff < 0 else "+"

        # Variable und Exponent
        if exp == 0:
            # Konstante
            var_str = ""
        elif exp == 1:
            # Linearer Term
            var_str = var
        else:
            # Höhere Potenzen
            var_str = f"{var}^{{{exp}}}"

        # Koeffizient formatieren
        if var_str == "":
            # Konstanter Term: immer Wert anzeigen
            if isinstance(abs_coeff, float) and abs_coeff != int(abs_coeff):
                coeff_str = f"{abs_coeff}"
            else:
                coeff_str = str(int(abs_coeff)) if isinstance(abs_coeff, float) else str(abs_coeff)
        elif abs_coeff == 1:
            # Koeffizient 1: weglassen
            coeff_str = ""
        else:
            if isinstance(abs_coeff, float) and abs_coeff != int(abs_coeff):
                coeff_str = f"{abs_coeff}"
            else:
                coeff_str = str(int(abs_coeff)) if isinstance(abs_coeff, float) else str(abs_coeff)

        terms.append((sign, coeff_str, var_str))

    if not terms:
        return "0"

    # Ersten Term zusammensetzen (kein "+" am Anfang)
    first_sign, first_coeff, first_var = terms[0]
    if first_sign == "-":
        result = "-" + first_coeff + first_var
    else:
        result = first_coeff + first_var

    # Weitere Terme
    for sign, coeff_str, var_str in terms[1:]:
        if sign == "-":
            result += " - " + coeff_str + var_str
        else:
            result += " + " + coeff_str + var_str

    return result


def matrix_to_latex(matrix_data: list, env: str = 'pmatrix') -> str:
    """
    @brief Konvertiert eine 2D-Liste in eine LaTeX-Matrix.
    @description
        Unterstützte Umgebungen:
        - 'pmatrix': Runde Klammern ( )
        - 'bmatrix': Eckige Klammern [ ]
        - 'vmatrix': Betragsstriche | |
        - 'matrix': Keine Klammern
        - 'Bmatrix': Geschweifte Klammern { }

        Format:
        \begin{pmatrix}
        a_{11} & a_{12} \\
        a_{21} & a_{22}
        \end{pmatrix}
    @param matrix_data 2D-Liste mit den Matrixeinträgen.
    @param env LaTeX-Umgebungsname (Standard: 'pmatrix').
    @return LaTeX-String der Matrix.
    @date 2026-03-09
    """
    lines = [r"\begin{" + env + "}"]

    for row in matrix_data:
        # Zeile: Elemente durch " & " trennen, Zeilenende " \\"
        row_str = " & ".join(number_to_latex(elem) if isinstance(elem, (int, float, complex))
                             else str(elem)
                             for elem in row)
        lines.append(row_str + r" \\")

    # Letzte "\\" entfernen (optionales LaTeX-Stilmittel)
    if len(lines) > 1:
        lines[-1] = lines[-1][:-3]  # " \\" entfernen

    lines.append(r"\end{" + env + "}")
    return "\n".join(lines)


def vector_to_latex(components: list, column: bool = True) -> str:
    """
    @brief Konvertiert eine Komponentenliste in einen LaTeX-Vektor.
    @description
        Spaltenvektor (column=True):
            \begin{pmatrix} a \\ b \\ c \end{pmatrix}
        Zeilenvektor (column=False):
            \begin{pmatrix} a & b & c \end{pmatrix}
    @param components Liste der Vektorkomponenten.
    @param column True für Spaltenvektor, False für Zeilenvektor.
    @return LaTeX-String des Vektors.
    @date 2026-03-09
    """
    # Elemente als LaTeX-Strings
    elem_strs = [number_to_latex(c) if isinstance(c, (int, float, complex))
                 else str(c)
                 for c in components]

    if column:
        # Spaltenvektor: eine Komponente pro Zeile
        inner = r" \\ ".join(elem_strs)
        return r"\begin{pmatrix} " + inner + r" \end{pmatrix}"
    else:
        # Zeilenvektor: alle Komponenten nebeneinander
        inner = " & ".join(elem_strs)
        return r"\begin{pmatrix} " + inner + r" \end{pmatrix}"


def equation_to_latex(lhs: str, rhs: str, eq_type: str = '=') -> str:
    """
    @brief Erstellt einen LaTeX-Gleichungsstring.
    @description
        Unterstützte eq_type-Werte:
        - '='       → lhs = rhs
        - '<='      → lhs \leq rhs
        - '>='      → lhs \geq rhs
        - '!='      → lhs \neq rhs
        - '~='      → lhs \approx rhs
        - 'iff'     → lhs \iff rhs
        - 'implies' → lhs \implies rhs
        - '<'       → lhs < rhs
        - '>'       → lhs > rhs
    @param lhs Linke Seite der Gleichung.
    @param rhs Rechte Seite der Gleichung.
    @param eq_type Art der Beziehung.
    @return LaTeX-Gleichungsstring.
    @date 2026-03-09
    """
    # Mapping von eq_type zu LaTeX-Symbol
    symbols = {
        '=':        '=',
        '<=':       r'\leq',
        '>=':       r'\geq',
        '!=':       r'\neq',
        '~=':       r'\approx',
        'iff':      r'\iff',
        'implies':  r'\implies',
        '<':        '<',
        '>':        '>',
        'in':       r'\in',
        'subset':   r'\subseteq',
    }

    symbol = symbols.get(eq_type, '=')
    return f"{lhs} {symbol} {rhs}"


def sympy_to_latex(expr) -> str:
    """
    @brief Konvertiert einen SymPy-Ausdruck in LaTeX.
    @description
        Nutzt sp.latex() für die Konvertierung. Gibt den LaTeX-String
        zurück, der direkt in Matheumgebungen verwendet werden kann.
    @param expr SymPy-Ausdruck.
    @return LaTeX-String des Ausdrucks.
    @date 2026-03-09
    """
    return sp.latex(expr)


def solution_to_latex(result: dict, title: str = '') -> str:
    """
    @brief Konvertiert ein Ergebnis-Dictionary in einen LaTeX-Block.
    @description
        Erstellt eine align*-Umgebung mit Key = Value Gleichungen.
        Optionaler Titel wird als Kommentar vorangestellt.
    @param result Dictionary mit Schlüssel-Wert-Paaren der Ergebnisse.
    @param title Optionaler Titel für den Block.
    @return LaTeX-String mit align*-Umgebung.
    @date 2026-03-09
    """
    lines = []

    if title:
        lines.append(r"\textbf{" + title + r"}")
        lines.append("")

    lines.append(r"\begin{align*}")

    items = list(result.items())
    for i, (key, value) in enumerate(items):
        # Wert in LaTeX konvertieren
        if isinstance(value, (int, float, complex)):
            val_str = number_to_latex(value)
        elif isinstance(value, list):
            val_str = str(value)
        else:
            val_str = str(value)

        # Letzte Zeile ohne \\
        if i < len(items) - 1:
            lines.append(f"  {key} &= {val_str} \\\\")
        else:
            lines.append(f"  {key} &= {val_str}")

    lines.append(r"\end{align*}")
    return "\n".join(lines)


def integral_to_latex(integrand: str, var: str, lower=None, upper=None) -> str:
    """
    @brief Erstellt einen LaTeX-Integral-String.
    @description
        Bestimmtes Integral (lower und upper angegeben):
            \int_{lower}^{upper} integrand \, d{var}
        Unbestimmtes Integral (ohne Grenzen):
            \int integrand \, d{var}
    @param integrand Der Integrand als String.
    @param var Integrationsvariable.
    @param lower Untere Integrationsgrenze (optional).
    @param upper Obere Integrationsgrenze (optional).
    @return LaTeX-Integralstring.
    @date 2026-03-09
    """
    if lower is not None and upper is not None:
        # Bestimmtes Integral
        lower_str = str(lower)
        upper_str = str(upper)
        return r"\int_{" + lower_str + r"}^{" + upper_str + r"} " + integrand + r" \, d" + var
    else:
        # Unbestimmtes Integral
        return r"\int " + integrand + r" \, d" + var


def sum_to_latex(term: str, index: str, lower: str, upper: str) -> str:
    """
    @brief Erstellt einen LaTeX-Summen-String.
    @description
        Erzeugt: \sum_{index=lower}^{upper} term
    @param term Der Summand als String.
    @param index Laufindex-Variable.
    @param lower Untere Summationsgrenze.
    @param upper Obere Summationsgrenze.
    @return LaTeX-Summenstring.
    @date 2026-03-09
    """
    return r"\sum_{" + index + "=" + lower + r"}^{" + upper + r"} " + term


def limit_to_latex(expr: str, var: str, point: str, direction: str = '') -> str:
    """
    @brief Erstellt einen LaTeX-Grenzwert-String.
    @description
        Ohne Richtung:  \lim_{var \to point} expr
        Mit '+':        \lim_{var \to point^+} expr
        Mit '-':        \lim_{var \to point^-} expr
    @param expr Der Ausdruck dessen Grenzwert berechnet wird.
    @param var Grenzwertvariable.
    @param point Punkt gegen den die Variable strebt.
    @param direction Richtung: '+', '-' oder '' für beidseitig.
    @return LaTeX-Grenzwertstring.
    @date 2026-03-09
    """
    if direction in ('+', '-'):
        sub = r"\to " + str(point) + "^{" + direction + "}"
    else:
        sub = r"\to " + str(point)

    return r"\lim_{" + var + " " + sub + r"} " + expr


def derivative_to_latex(func: str, var: str, order: int = 1) -> str:
    """
    @brief Erstellt einen LaTeX-Ableitungsstring.
    @description
        Erste Ableitung (order=1):  \frac{d}{dx} f(x)
        n-te Ableitung (order>1):   \frac{d^n}{dx^n} f(x)
    @param func Funktion als String.
    @param var Ableitungsvariable.
    @param order Ordnung der Ableitung.
    @return LaTeX-Ableitungsstring.
    @date 2026-03-09
    """
    if order == 1:
        return r"\frac{d}{d" + var + r"} " + func
    else:
        return r"\frac{d^{" + str(order) + r"}}{d" + var + r"^{" + str(order) + r"}} " + func


def theorem_to_latex(name: str, statement: str, proof: str = '') -> str:
    """
    @brief Erstellt einen LaTeX-Theorem-Block.
    @description
        Format:
        \begin{theorem}[name]
        statement
        \end{theorem}
        \begin{proof}
        proof
        \end{proof}

        Der proof-Block wird nur erzeugt wenn proof nicht leer ist.
        Benötigt das amsthm-Paket in der LaTeX-Präambel.
    @param name Name oder Bezeichnung des Theorems.
    @param statement Aussage/Formulierung des Theorems.
    @param proof Beweis des Theorems (optional).
    @return LaTeX-String mit theorem (und ggf. proof) Umgebung.
    @date 2026-03-09
    """
    lines = []
    lines.append(r"\begin{theorem}[" + name + r"]")
    lines.append(statement)
    lines.append(r"\end{theorem}")

    if proof:
        lines.append(r"\begin{proof}")
        lines.append(proof)
        lines.append(r"\end{proof}")

    return "\n".join(lines)


def full_document(content: str,
                  title: str = 'Mathematische Ergebnisse',
                  author: str = 'Michael Fuhrmann',
                  packages: list = None) -> str:
    """
    @brief Erstellt ein vollständiges LaTeX-Dokument.
    @description
        Erzeugt ein gültiges LaTeX-Dokument mit:
        - documentclass{article}
        - Standard mathematischen Paketen (amsmath, amsfonts, amssymb, amsthm)
        - Optionalen zusätzlichen Paketen
        - Titel und Autor
        - Übergebenem Inhalt

        Das Dokument kann direkt mit pdflatex, lualatex oder xelatex kompiliert werden.
    @param content LaTeX-Inhalt des Dokuments (zwischen \begin{document} und \end{document}).
    @param title Titel des Dokuments.
    @param author Autor des Dokuments.
    @param packages Zusätzliche Pakete als Liste von Strings, z.B. ['graphicx', 'hyperref'].
    @return Vollständiger LaTeX-Dokumentstring.
    @date 2026-03-09
    """
    # Standard-Pakete
    standard_packages = ['amsmath', 'amsfonts', 'amssymb', 'amsthm', 'geometry']
    all_packages = standard_packages.copy()

    if packages:
        for pkg in packages:
            if pkg not in all_packages:
                all_packages.append(pkg)

    # Pakete als LaTeX-Zeilen
    package_lines = "\n".join(r"\usepackage{" + pkg + r"}" for pkg in all_packages)

    doc = r"""\documentclass{article}
""" + package_lines + r"""
\geometry{a4paper, margin=2.5cm}

\title{""" + title + r"""}
\author{""" + author + r"""}
\date{\today}

\begin{document}
\maketitle

""" + content + r"""

\end{document}"""

    return doc


def save_latex(content: str, filepath: str) -> None:
    """
    @brief Speichert einen LaTeX-String in eine .tex-Datei.
    @description
        Schreibt den übergebenen LaTeX-Inhalt in eine Datei.
        Die Datei wird mit UTF-8-Kodierung gespeichert.
        Falls die Datei bereits existiert, wird sie überschrieben.
    @param content LaTeX-Inhalt als String.
    @param filepath Pfad zur Zieldatei (empfohlen mit .tex Endung).
    @return None
    @date 2026-03-09
    """
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
