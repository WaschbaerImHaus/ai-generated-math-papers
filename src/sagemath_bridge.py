"""
@file sagemath_bridge.py
@brief SageMath-Interoperabilitäts-Brücke.
@description
    Stellt Funktionen für die Konvertierung zwischen SymPy und SageMath
    bereit. Alle Funktionen arbeiten auch ohne installiertes SageMath
    (Fallback-Modus).

    Wenn SageMath installiert ist:
    - to_sage(): Konvertiert SymPy-Ausdrücke zu SageMath
    - from_sage(): Konvertiert SageMath-Ausdrücke zurück zu SymPy

    Ohne SageMath:
    - Alle Konvertierungen geben SymPy-Ausdrücke zurück
    - export_to_sage_format() erstellt SageMath-kompatible Python-Skripte

@author Michael Fuhrmann
@date 2026-03-11
@lastModified 2026-03-11
"""

import sympy as sp
from typing import Any, Optional


def is_sage_available() -> bool:
    """
    @brief Prüft ob SageMath installiert und importierbar ist.
    @description
        Versucht das sage-Modul zu importieren. Gibt True zurück
        wenn SageMath verfügbar ist, sonst False.

        SageMath ist ein Open-Source-Mathematiksystem, das viele
        Bibliotheken (SymPy, NumPy, FLINT, etc.) integriert.

    @return: True wenn SageMath verfügbar, sonst False.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    try:
        import sage.all  # type: ignore[import]
        return True
    except ImportError:
        return False


def to_sage(sympy_expr: Any) -> Any:
    """
    @brief Konvertiert einen SymPy-Ausdruck zu SageMath.
    @description
        Wenn SageMath verfügbar ist, wird der Ausdruck konvertiert.
        SageMath kann den Ausdruck dann mit seinen eigenen Methoden
        weiterverarbeiten (z.B. prove(), plot(), etc.).

        Wenn SageMath nicht verfügbar ist, wird der SymPy-Ausdruck
        unverändert zurückgegeben (Fallback-Modus).

        Konvertierungsweg:
        SymPy → LaTeX-String → SageMath (via sage_eval oder SR())

    @param sympy_expr: SymPy-Ausdruck oder konvertierbarer Wert.
    @return: SageMath-Ausdruck oder SymPy-Fallback.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Eingabe zu SymPy konvertieren falls nötig
    if not isinstance(sympy_expr, sp.Basic):
        try:
            sympy_expr = sp.sympify(sympy_expr)
        except Exception:
            return sympy_expr

    # Versuche SageMath zu nutzen
    if is_sage_available():
        try:
            import sage.all as sage  # type: ignore[import]
            # Konvertierung via SymPy-zu-Sage-Schnittstelle
            return sage.SR(sympy_expr._sage_())
        except Exception:
            pass  # Fallback bei Konvertierungsfehler

    # Fallback: SymPy-Ausdruck zurückgeben
    return sympy_expr


def from_sage(sage_expr: Any) -> sp.Expr:
    """
    @brief Konvertiert einen SageMath-Ausdruck zurück zu SymPy.
    @description
        Wenn SageMath verfügbar ist, wird der Ausdruck nach SymPy
        konvertiert. Andernfalls wird der Ausdruck unverändert
        zurückgegeben (wird als SymPy-kompatibel angenommen).

    @param sage_expr: SageMath-Ausdruck.
    @return: SymPy-Ausdruck.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    if is_sage_available():
        try:
            # SageMath-Ausdrücke haben eine _sympy_()-Methode
            if hasattr(sage_expr, '_sympy_'):
                return sage_expr._sympy_()
        except Exception:
            pass

    # Fallback: sympify (funktioniert bei String-Darstellungen)
    try:
        return sp.sympify(sage_expr)
    except Exception:
        return sage_expr


def export_to_sage_format(data_dict: dict, filename: Optional[str] = None) -> str:
    """
    @brief Exportiert Ergebnisse als SageMath-kompatibles Python-Skript.
    @description
        Erstellt ein Python-Skript, das sowohl mit SageMath als auch
        mit reinem SymPy ausgeführt werden kann.

        Struktur des erzeugten Skripts:
        1. Imports (SageMath-Kompatibilität prüfen)
        2. Symbolische Variablen definieren
        3. Alle Ergebnisse als SymPy-Ausdrücke
        4. SageMath-spezifische Operationen (wenn verfügbar)

        Das Skript funktioniert AUCH ohne installiertes SageMath,
        da es auf SymPy als Fallback zurückgreift.

    @param data_dict: Dictionary mit Ergebnis-Namen und SymPy-Ausdrücken.
    @param filename: Optionaler Dateipfad zum Speichern (None = kein Speichern).
    @return: Python-Skript als String.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Skript-Kopf: SageMath-Kompatibilitäts-Prüfung
    lines = [
        "#!/usr/bin/env python3",
        '"""',
        "SageMath-kompatibles Exportskript – erzeugt von specialist-maths",
        "Läuft mit SageMath (sage script.py) oder Python3 (python3 script.py)",
        '"""',
        "",
        "# Versuche SageMath zu importieren; Fallback auf SymPy",
        "try:",
        "    from sage.all import *",
        "    SAGE_AVAILABLE = True",
        "except ImportError:",
        "    import sympy as sp",
        "    from sympy import *",
        "    SAGE_AVAILABLE = False",
        "",
        "print(f'SageMath verfügbar: {SAGE_AVAILABLE}')",
        "",
        "# Symbolische Variablen",
        "if SAGE_AVAILABLE:",
        "    x, y, z, t, n, k = var('x y z t n k')",
        "else:",
        "    x, y, z, t, n, k = sp.symbols('x y z t n k')",
        "",
        "# --- Exportierte Ergebnisse ---",
    ]

    # Jedes Ergebnis ins Skript einfügen
    for name, value in data_dict.items():
        # Sicheren Python-Variablennamen erzeugen
        safe_name = _make_safe_varname(name)

        try:
            # SymPy-Ausdruck zu String konvertieren
            if isinstance(value, sp.Basic):
                expr_str = str(value)
                latex_str = sp.latex(value)
            elif isinstance(value, (int, float)):
                expr_str = str(value)
                latex_str = str(value)
            elif isinstance(value, (list, tuple)):
                expr_str = str(value)
                latex_str = str(value)
            else:
                expr_str = repr(value)
                latex_str = str(value)

            lines.append(f"")
            lines.append(f"# Ergebnis: {name}")
            lines.append(f"# LaTeX: {latex_str}")
            lines.append(f"{safe_name} = sympify('{expr_str}') if not SAGE_AVAILABLE "
                         f"else SR('{expr_str}')")
            lines.append(f"print(f'{name} = {{{safe_name}}}')")
        except Exception as e:
            lines.append(f"# Fehler beim Export von '{name}': {e}")
            lines.append(f"{safe_name} = None")

    # Abschluss
    lines.extend([
        "",
        "# SageMath-spezifische Ausgaben (nur wenn SageMath verfügbar)",
        "if SAGE_AVAILABLE:",
        "    print('SageMath-Modus aktiv: Erweiterte Operationen verfügbar')",
        "    # Beispiel: show(x, viewer='ascii')",
    ])

    # Skript zusammenfügen
    script = "\n".join(lines)

    # Optional in Datei speichern
    if filename:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(script)

    return script


def export_to_sage_string(sympy_expr: Any) -> str:
    """
    @brief Konvertiert einen SymPy-Ausdruck zu einem SageMath-kompatiblen String.
    @description
        Erzeugt eine String-Darstellung des SymPy-Ausdrucks, die direkt in
        SageMath eingegeben werden kann (oder in einem Python-Skript, das
        SageMath-Syntax verwendet).

        Konvertierungsregeln:
        - SymPy-Ausdrücke → str(expr) (kompatibel mit SageMath)
        - Spezielle Konstanten: sympy.pi → „pi", sympy.E → „e"
        - Funktionen: sympy.sin → „sin", sympy.cos → „cos" usw.
        - Unbekannte Ausdrücke → repr()

        Diese Funktion funktioniert OHNE installiertes SageMath.
        Sie erstellt nur den String; die tatsächliche Auswertung
        in SageMath muss separat erfolgen.

    @param sympy_expr: SymPy-Ausdruck (sp.Basic) oder konvertierbarer Wert.
    @return: SageMath-kompatibler String.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    # Eingabe zu SymPy konvertieren falls nötig
    if not isinstance(sympy_expr, sp.Basic):
        try:
            sympy_expr = sp.sympify(sympy_expr)
        except Exception:
            # Fallback: direkte String-Konvertierung
            return str(sympy_expr)

    # SageMath-spezifische Ausgabe via SymPy-Printer
    # SymPy's str()-Repräsentation ist weitgehend SageMath-kompatibel
    try:
        # Versuche SageMath-spezifische Darstellung via sage_repr
        if is_sage_available():
            sage_expr = to_sage(sympy_expr)
            return str(sage_expr)
        else:
            # Ohne SageMath: SymPy-String (weitgehend kompatibel)
            return str(sympy_expr)
    except Exception:
        # Sicherer Fallback
        return str(sympy_expr)


def import_from_sage_string(sage_str: str) -> sp.Expr:
    """
    @brief Parst einen SageMath-Output-String und gibt einen SymPy-Ausdruck zurück.
    @description
        Versucht einen String, der SageMath-Syntax enthält, in einen
        SymPy-Ausdruck zu konvertieren.

        Strategie:
        1. Falls SageMath verfügbar: sage_eval() → _sympy_()
        2. Fallback: sympy.sympify() mit angepasster Lokalmenge
           (Warnung wird geloggt)
        3. Letzter Fallback: sympy.Symbol(sage_str) wenn alles fehlschlägt

        Unterstützte SageMath-Syntax:
        - Variablen: x, y, z, t (werden als SymPy-Symbols interpretiert)
        - Operatoren: +, -, *, /, ^, ** (^ wird zu ** konvertiert)
        - Funktionen: sin, cos, tan, exp, log, sqrt
        - Konstanten: pi, e, I (werden zu sp.pi, sp.E, sp.I)

    @param sage_str: SageMath-String (z.B. "x^2 + sin(x)").
    @return: SymPy-Ausdruck.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    import warnings

    # Leerer String → symbolischer Ausdruck 0
    if not sage_str or not sage_str.strip():
        return sp.Integer(0)

    # Strategie 1: SageMath direkt nutzen (wenn verfügbar)
    if is_sage_available():
        try:
            import sage.all as sage  # type: ignore[import]
            sage_expr = sage.sage_eval(sage_str)
            # Zu SymPy konvertieren
            return from_sage(sage_expr)
        except Exception:
            pass  # Fallback auf SymPy

    # Warnung ausgeben: Kein SageMath verfügbar, nutze SymPy-Fallback
    warnings.warn(
        "SageMath nicht verfügbar – parse_from_sage_string() nutzt SymPy-Parser "
        "als Fallback. Komplexe SageMath-Syntax könnte nicht korrekt interpretiert werden.",
        UserWarning,
        stacklevel=2
    )

    # SageMath → Python-Konvertierungen
    # ^ ist Potenz in SageMath, ** in Python/SymPy
    python_str = sage_str.replace('^', '**')

    # Lokale Symbole und Funktionen für sympify bereitstellen
    local_dict = {
        'x': sp.Symbol('x'),
        'y': sp.Symbol('y'),
        'z': sp.Symbol('z'),
        't': sp.Symbol('t'),
        'n': sp.Symbol('n'),
        'k': sp.Symbol('k'),
        'pi': sp.pi,
        'e':  sp.E,
        'I':  sp.I,
        'sin': sp.sin,
        'cos': sp.cos,
        'tan': sp.tan,
        'exp': sp.exp,
        'log': sp.log,
        'sqrt': sp.sqrt,
        'abs': sp.Abs,
        'conjugate': sp.conjugate,
    }

    # Strategie 2: SymPy-Parser mit angepasster Lokalmenge
    try:
        return sp.sympify(python_str, locals=local_dict)
    except Exception:
        pass

    # Strategie 3: Direkte sympify ohne Hilfsvariablen
    try:
        return sp.sympify(python_str)
    except Exception:
        pass

    # Letzter Fallback: Als Symbol-Name behandeln
    return sp.Symbol(sage_str.strip())


def _make_safe_varname(name: str) -> str:
    """
    @brief Konvertiert einen String in einen sicheren Python-Variablennamen.
    @description
        Ersetzt Sonderzeichen durch Unterstriche und stellt sicher,
        dass der Name mit einem Buchstaben oder Unterstrich beginnt.

    @param name: Zu konvertierender Name.
    @return: Sicherer Python-Variablenname.
    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """
    import re
    # Alle Nicht-Alphanumerischen-Zeichen durch Unterstriche ersetzen
    safe = re.sub(r'[^a-zA-Z0-9_]', '_', name)

    # Name darf nicht mit Ziffer beginnen
    if safe and safe[0].isdigit():
        safe = 'var_' + safe

    # Leerer Name → Standardname
    if not safe:
        safe = 'result'

    return safe
