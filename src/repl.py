"""
@file repl.py
@brief Interaktiver REPL-Modus für mathematische Berechnungen.
@description
    Stellt eine Read-Eval-Print-Loop (REPL) Schnittstelle bereit,
    über die alle mathematischen Module interaktiv genutzt werden können.
    Unterstützt auch Jupyter-Notebook-Kompatibilität.

    Verfügbare Befehle:
    - prime N        : Primzahltest für N
    - factor N       : Primfaktorzerlegung von N
    - gcd A B        : Größter gemeinsamer Teiler
    - lcm A B        : Kleinstes gemeinsames Vielfaches
    - solve A B C    : Quadratische Gleichung ax²+bx+c=0
    - derive EXPR at X : Numerische Ableitung
    - integrate F A B  : Numerisches Integral
    - limit EXPR VAR PT: Symbolischer Grenzwert
    - eigenvalues MAT  : Eigenwerte einer Matrix
    - fft LIST         : Fast Fourier Transformation
    - goldbach N       : Goldbach-Zerlegung
    - collatz N        : Collatz-Stoppzeit
    - help             : Alle Befehle anzeigen
    - exit/quit        : REPL beenden

@author Kurt Ingwer
@date 2026-03-09
@version 1.0.0
"""

import math
import json
import sys
import os
from typing import Any, Optional

# Pfad zum src-Verzeichnis hinzufügen
_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

# Mathematische Module importieren
import algebra
import analysis
import linear_algebra
import fourier
import proof_theory
import statistics_math


def is_jupyter() -> bool:
    """
    @brief Prüft ob das Programm in einem Jupyter-Notebook läuft.
    @description
        Versucht get_ipython() aufzurufen. Falls es existiert und eine
        Jupyter-Umgebung zurückgibt (ZMQInteractiveShell), läuft das
        Programm in einem Jupyter-Notebook.
    @return True wenn in Jupyter, sonst False.
    @date 2026-03-09
    """
    try:
        # get_ipython ist nur in IPython/Jupyter verfügbar
        shell = get_ipython().__class__.__name__  # type: ignore[name-defined]
        if shell == 'ZMQInteractiveShell':
            # Jupyter Notebook oder JupyterLab
            return True
        elif shell == 'TerminalInteractiveShell':
            # IPython Terminal
            return False
        else:
            return False
    except NameError:
        # Normales Python-Skript
        return False


def create_jupyter_notebook(output_path: str) -> None:
    """
    @brief Erstellt eine Jupyter-Notebook-Datei (.ipynb) mit Demo-Zellen.
    @description
        Erzeugt ein Jupyter-Notebook im nbformat 4.4 JSON-Format,
        das Demo-Zellen für alle verfügbaren mathematischen Module enthält.
        Falls nbformat installiert ist, wird es genutzt; andernfalls wird
        das JSON direkt geschrieben.
    @param output_path Pfad zur .ipynb-Ausgabedatei.
    @return None
    @date 2026-03-09
    """
    # Notebook-Struktur als Dictionary
    notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
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
        "cells": []
    }

    def make_code_cell(source: str) -> dict:
        """Erstellt eine Code-Zelle."""
        return {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": source
        }

    def make_markdown_cell(source: str) -> dict:
        """Erstellt eine Markdown-Zelle."""
        return {
            "cell_type": "markdown",
            "metadata": {},
            "source": source
        }

    # Zellen hinzufügen
    cells = notebook["cells"]

    # Titel
    cells.append(make_markdown_cell(
        "# specialist-maths - Demo Notebook\n\n"
        "Interaktive Demonstration aller mathematischen Module."
    ))

    # Setup
    cells.append(make_code_cell(
        "import sys, os\n"
        "sys.path.insert(0, '../src')\n\n"
        "import algebra\n"
        "import analysis\n"
        "import linear_algebra\n"
        "import fourier\n"
        "import proof_theory\n"
        "import statistics_math\n"
        "import latex_export as lx\n"
        "print('Module geladen')"
    ))

    # Algebra
    cells.append(make_markdown_cell("## Algebra"))
    cells.append(make_code_cell(
        "# Primzahltest\n"
        "print(algebra.is_prime(97))\n\n"
        "# Primfaktorzerlegung\n"
        "print(algebra.prime_factorization(360))\n\n"
        "# GGT und KGV\n"
        "print(algebra.gcd(48, 18))\n"
        "print(algebra.lcm(12, 18))\n\n"
        "# Quadratische Gleichung\n"
        "print(algebra.solve_quadratic(1, -5, 6))"
    ))

    # Analysis
    cells.append(make_markdown_cell("## Analysis"))
    cells.append(make_code_cell(
        "import math\n\n"
        "# Numerische Ableitung von sin(x) bei x=0\n"
        "deriv = analysis.numerical_derivative(math.sin, 0)\n"
        "print(f'sin\\'(0) = {deriv:.6f}')\n\n"
        "# Numerisches Integral\n"
        "result = analysis.numerical_integral(math.sin, 0, math.pi)\n"
        "print(f'∫sin(x)dx von 0 bis π = {result:.6f}')\n\n"
        "# Symbolischer Grenzwert\n"
        "lim = analysis.symbolic_limit('sin(x)/x', 'x', 0)\n"
        "print(f'lim(sin(x)/x, x→0) = {lim}')"
    ))

    # Lineare Algebra
    cells.append(make_markdown_cell("## Lineare Algebra"))
    cells.append(make_code_cell(
        "# Vektor-Operationen\n"
        "v1 = linear_algebra.Vector([1, 2, 3])\n"
        "v2 = linear_algebra.Vector([4, 5, 6])\n"
        "print(f'Skalarprodukt: {v1.dot(v2)}')\n"
        "print(f'Kreuzprodukt: {v1.cross(v2).components}')\n\n"
        "# Matrix-Eigenwerte\n"
        "M = linear_algebra.Matrix([[4, 1], [2, 3]])\n"
        "print(f'Eigenwerte: {M.eigenvalues()}')"
    ))

    # FFT
    cells.append(make_markdown_cell("## Fourier-Transformation"))
    cells.append(make_code_cell(
        "# FFT eines Signals\n"
        "signal = [1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.0]\n"
        "fft_result = fourier.fft(signal)\n"
        "magnitudes = [abs(c) for c in fft_result]\n"
        "print('FFT-Amplituden:', [f'{m:.3f}' for m in magnitudes])"
    ))

    # LaTeX-Export
    cells.append(make_markdown_cell("## LaTeX-Export"))
    cells.append(make_code_cell(
        "# Zahl in LaTeX\n"
        "print(lx.number_to_latex(3.14159))\n"
        "print(lx.number_to_latex(math.pi))\n\n"
        "# Polynom in LaTeX\n"
        "print(lx.polynomial_to_latex([1, -3, 2]))\n\n"
        "# Matrix in LaTeX\n"
        "print(lx.matrix_to_latex([[1, 2], [3, 4]]))"
    ))

    # Notebook speichern
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=1, ensure_ascii=False)


# =============================================================================
# KLASSE: MathREPL
# =============================================================================

class MathREPL:
    """
    @brief Interaktiver REPL (Read-Eval-Print-Loop) für Mathematik.
    @description
        Bietet eine einfache Kommandozeilen-Schnittstelle zu allen
        mathematischen Modulen. Befehle werden als Strings eingegeben
        und das Ergebnis als String zurückgegeben.

    @example
        repl = MathREPL()
        repl.run()  # Startet interaktive Schleife
    @date 2026-03-09
    """

    def __init__(self):
        """
        @brief Initialisiert den REPL und lädt alle mathematischen Module.
        @description
            Stellt Referenzen auf alle Module bereit und definiert
            die Befehlstabelle für evaluate().
        @date 2026-03-09
        """
        # Module als Attribute speichern
        self._algebra = algebra
        self._analysis = analysis
        self._linalg = linear_algebra
        self._fourier = fourier
        self._proof = proof_theory
        self._stats = statistics_math

        # Befehlstabelle: Name → (Handler-Methode, Beschreibung)
        self._commands = self._build_command_table()

    def _build_command_table(self) -> dict:
        """
        @brief Erstellt die Tabelle der verfügbaren Befehle.
        @return Dictionary {Befehlsname: Beschreibung}.
        @date 2026-03-09
        """
        return {
            'prime':        'prime N          — Primzahltest für N',
            'factor':       'factor N         — Primfaktorzerlegung von N',
            'gcd':          'gcd A B          — Größter gemeinsamer Teiler von A und B',
            'lcm':          'lcm A B          — Kleinstes gemeinsames Vielfaches von A und B',
            'solve':        'solve A B C      — Quadratische Gleichung Ax²+Bx+C=0',
            'derive':       'derive FUNC at X — Numerische Ableitung (FUNC: sin/cos/exp/sq)',
            'integrate':    'integrate F A B  — Numerisches Integral von F von A bis B',
            'limit':        'limit EXPR VAR PT— Symbolischer Grenzwert',
            'eigenvalues':  'eigenvalues JSON — Eigenwerte einer Matrix (JSON-Format)',
            'fft':          'fft LIST         — Fast Fourier Transformation (JSON-Liste)',
            'goldbach':     'goldbach N       — Goldbach-Zerlegung von N (gerade N > 2)',
            'collatz':      'collatz N        — Collatz-Stoppzeit für N',
            'phi':          'phi N            — Euler-Phi-Funktion φ(N)',
            'modinv':       'modinv A M       — Modulares Inverses von A mod M',
            'sqrt':         'sqrt N           — Quadratwurzel von N',
            'help':         'help [BEFEHL]    — Alle Befehle oder Hilfe zu BEFEHL',
            'exit':         'exit             — REPL beenden',
            'quit':         'quit             — REPL beenden',
        }

    def _parse_function(self, func_name: str):
        """
        @brief Konvertiert einen Funktionsnamen in ein aufrufbares Objekt.
        @description
            Unterstützte Funktionsnamen: sin, cos, tan, exp, log, sqrt, sq (x²), id (x)
        @param func_name Name der Funktion als String.
        @return Aufrufbares Python-Objekt.
        @date 2026-03-09
        """
        func_map = {
            'sin':   math.sin,
            'cos':   math.cos,
            'tan':   math.tan,
            'exp':   math.exp,
            'log':   math.log,
            'sqrt':  math.sqrt,
            'sq':    lambda x: x * x,
            'id':    lambda x: x,
            'abs':   abs,
        }
        if func_name not in func_map:
            raise ValueError(f"Unbekannte Funktion '{func_name}'. "
                             f"Verfügbar: {', '.join(func_map.keys())}")
        return func_map[func_name]

    def evaluate(self, command: str) -> str:
        """
        @brief Wertet einen REPL-Befehl aus und gibt das Ergebnis zurück.
        @description
            Parst den Befehlsstring, ruft die entsprechende mathematische
            Funktion auf und gibt das Ergebnis als formatierten String zurück.
            Bei Fehlern wird eine informative Fehlermeldung zurückgegeben.
        @param command Befehlsstring (z.B. "prime 97", "gcd 48 18").
        @return Ergebnis als String.
        @date 2026-03-09
        """
        # Leerzeichen trimmen und leere Befehle abfangen
        command = command.strip()
        if not command:
            return ""

        # Befehl aufteilen
        parts = command.split()
        cmd = parts[0].lower()
        args = parts[1:]

        try:
            # -------------------------
            # Befehl: prime N
            # -------------------------
            if cmd == 'prime':
                if len(args) != 1:
                    return "Fehler: prime erwartet genau 1 Argument (Zahl)"
                n = int(args[0])
                result = self._algebra.is_prime(n)
                return f"is_prime({n}) = {result}"

            # -------------------------
            # Befehl: factor N
            # -------------------------
            elif cmd == 'factor':
                if len(args) != 1:
                    return "Fehler: factor erwartet genau 1 Argument"
                n = int(args[0])
                result = self._algebra.prime_factorization(n)
                # Lesbare Darstellung
                factors_str = " × ".join(
                    f"{p}^{e}" if e > 1 else str(p)
                    for p, e in sorted(result.items())
                )
                return f"prime_factorization({n}) = {result}\n  = {factors_str}"

            # -------------------------
            # Befehl: gcd A B
            # -------------------------
            elif cmd == 'gcd':
                if len(args) != 2:
                    return "Fehler: gcd erwartet 2 Argumente"
                a, b = int(args[0]), int(args[1])
                result = self._algebra.gcd(a, b)
                return f"gcd({a}, {b}) = {result}"

            # -------------------------
            # Befehl: lcm A B
            # -------------------------
            elif cmd == 'lcm':
                if len(args) != 2:
                    return "Fehler: lcm erwartet 2 Argumente"
                a, b = int(args[0]), int(args[1])
                result = self._algebra.lcm(a, b)
                return f"lcm({a}, {b}) = {result}"

            # -------------------------
            # Befehl: solve A B C
            # -------------------------
            elif cmd == 'solve':
                if len(args) != 3:
                    return "Fehler: solve erwartet 3 Koeffizienten (a b c)"
                a, b, c = float(args[0]), float(args[1]), float(args[2])
                roots = self._algebra.solve_quadratic(a, b, c)
                return f"solve_quadratic({a}x² + {b}x + {c} = 0) = {roots}"

            # -------------------------
            # Befehl: derive FUNC at X
            # -------------------------
            elif cmd == 'derive':
                # Format: "derive FUNC at X"
                if len(args) < 3 or args[1].lower() != 'at':
                    return "Fehler: Format: derive FUNC at X"
                func = self._parse_function(args[0])
                x = float(args[2])
                result = self._analysis.numerical_derivative(func, x)
                return f"d/dx {args[0]}(x) bei x={x} ≈ {result:.8f}"

            # -------------------------
            # Befehl: integrate F A B
            # -------------------------
            elif cmd == 'integrate':
                if len(args) != 3:
                    return "Fehler: integrate erwartet 3 Argumente: FUNC A B"
                func = self._parse_function(args[0])
                a, b = float(args[1]), float(args[2])
                result = self._analysis.numerical_integral(func, a, b)
                return f"∫{args[0]}(x)dx von {a} bis {b} ≈ {result:.8f}"

            # -------------------------
            # Befehl: limit EXPR VAR PT
            # -------------------------
            elif cmd == 'limit':
                if len(args) < 3:
                    return "Fehler: limit erwartet 3 Argumente: EXPR VAR PUNKT"
                expr_str = args[0]
                var = args[1]
                # Punkt kann "inf", "pi", eine Zahl sein
                pt_str = args[2]
                try:
                    point = float(pt_str)
                except ValueError:
                    point = pt_str  # Symbolisch übergeben (z.B. "oo", "pi")
                result = self._analysis.symbolic_limit(expr_str, var, point)
                return f"lim({expr_str}, {var} → {pt_str}) = {result}"

            # -------------------------
            # Befehl: eigenvalues JSON-Matrix
            # -------------------------
            elif cmd == 'eigenvalues':
                if not args:
                    return "Fehler: eigenvalues erwartet eine Matrix als JSON"
                # Alle Args als JSON zusammensetzen
                json_str = " ".join(args)
                try:
                    matrix_data = json.loads(json_str)
                except json.JSONDecodeError:
                    return f"Fehler: Ungültiges JSON-Format für Matrix: {json_str}"
                M = self._linalg.Matrix(matrix_data)
                eigenvals = M.eigenvalues()
                return f"Eigenwerte = {[f'{v:.6f}' if isinstance(v, float) else str(v) for v in eigenvals]}"

            # -------------------------
            # Befehl: fft JSON-Liste
            # -------------------------
            elif cmd == 'fft':
                if not args:
                    return "Fehler: fft erwartet eine Liste als JSON"
                json_str = " ".join(args)
                try:
                    signal = json.loads(json_str)
                except json.JSONDecodeError:
                    return f"Fehler: Ungültiges JSON-Format: {json_str}"
                result = self._fourier.fft(signal)
                magnitudes = [abs(c) for c in result]
                return (f"FFT({signal}) =\n"
                        f"  Amplituden: {[f'{m:.4f}' for m in magnitudes]}")

            # -------------------------
            # Befehl: goldbach N
            # -------------------------
            elif cmd == 'goldbach':
                if len(args) != 1:
                    return "Fehler: goldbach erwartet 1 Argument"
                n = int(args[0])
                result = self._proof.goldbach_decomposition(n)
                if result:
                    return f"goldbach({n}) = {result[0]} + {result[1]}"
                else:
                    return f"goldbach({n}) = Keine Zerlegung gefunden"

            # -------------------------
            # Befehl: collatz N
            # -------------------------
            elif cmd == 'collatz':
                if len(args) != 1:
                    return "Fehler: collatz erwartet 1 Argument"
                n = int(args[0])
                result = self._proof.collatz_stopping_time(n)
                return f"collatz_stopping_time({n}) = {result}"

            # -------------------------
            # Befehl: phi N
            # -------------------------
            elif cmd == 'phi':
                if len(args) != 1:
                    return "Fehler: phi erwartet 1 Argument"
                n = int(args[0])
                result = self._algebra.euler_phi(n)
                return f"euler_phi({n}) = {result}"

            # -------------------------
            # Befehl: modinv A M
            # -------------------------
            elif cmd == 'modinv':
                if len(args) != 2:
                    return "Fehler: modinv erwartet 2 Argumente"
                a, m = int(args[0]), int(args[1])
                result = self._algebra.mod_inverse(a, m)
                return f"mod_inverse({a}, {m}) = {result}"

            # -------------------------
            # Befehl: sqrt N
            # -------------------------
            elif cmd == 'sqrt':
                if len(args) != 1:
                    return "Fehler: sqrt erwartet 1 Argument"
                n = float(args[0])
                result = math.sqrt(n)
                return f"sqrt({n}) = {result:.8f}"

            # -------------------------
            # Befehl: help
            # -------------------------
            elif cmd == 'help':
                if args:
                    # Hilfe zu einem bestimmten Befehl
                    target = args[0].lower()
                    if target in self._commands:
                        return f"Hilfe: {self._commands[target]}"
                    else:
                        return f"Unbekannter Befehl: {target}"
                else:
                    return self.help()

            # -------------------------
            # Befehl: exit/quit
            # -------------------------
            elif cmd in ('exit', 'quit'):
                return "__EXIT__"

            else:
                return (f"Unbekannter Befehl: '{cmd}'. "
                        f"Tippe 'help' für eine Befehlsübersicht.")

        except ValueError as e:
            return f"Eingabefehler: {e}"
        except Exception as e:
            return f"Fehler bei '{command}': {type(e).__name__}: {e}"

    def available_commands(self) -> dict:
        """
        @brief Gibt alle verfügbaren Befehle mit Beschreibungen zurück.
        @return Dictionary {Befehlsname: Beschreibungsstring}.
        @date 2026-03-09
        """
        return dict(self._commands)

    def help(self, command: str = '') -> str:
        """
        @brief Zeigt Hilfe zu einem Befehl oder allen Befehlen an.
        @description
            Ohne Argument: Listet alle verfügbaren Befehle auf.
            Mit Argument: Zeigt detaillierte Hilfe zum genannten Befehl.
        @param command Name des Befehls für detaillierte Hilfe (optional).
        @return Hilfetext als String.
        @date 2026-03-09
        """
        if command:
            cmd = command.lower()
            if cmd in self._commands:
                return f"Hilfe zu '{cmd}':\n  {self._commands[cmd]}"
            else:
                return f"Unbekannter Befehl: '{command}'"

        # Alle Befehle anzeigen
        lines = ["Verfügbare Befehle:", "=" * 50]
        for name, desc in self._commands.items():
            lines.append(f"  {desc}")
        lines.append("=" * 50)
        lines.append("Beispiele:")
        lines.append("  prime 97")
        lines.append("  gcd 48 18")
        lines.append("  solve 1 -5 6")
        lines.append("  derive sin at 1.5707963")
        lines.append("  integrate sin 0 3.14159")
        lines.append("  eigenvalues [[1,2],[3,4]]")
        lines.append("  fft [1,2,3,4,5,6,7,8]")
        return "\n".join(lines)

    def run(self):
        """
        @brief Startet die interaktive REPL-Schleife.
        @description
            Liest Befehle von der Standardeingabe, wertet sie aus und
            gibt das Ergebnis auf der Standardausgabe aus.
            Beenden mit 'exit', 'quit' oder Ctrl+C.
        @return None
        @date 2026-03-09
        """
        print("=" * 60)
        print("  specialist-maths REPL v1.0")
        print("  Tippe 'help' für Befehle, 'exit' zum Beenden")
        print("=" * 60)

        while True:
            try:
                # Eingabe lesen
                user_input = input("math> ").strip()

                if not user_input:
                    continue

                # Auswerten
                result = self.evaluate(user_input)

                # Beendigung prüfen
                if result == "__EXIT__":
                    print("Auf Wiedersehen!")
                    break

                # Ergebnis ausgeben
                print(result)

            except KeyboardInterrupt:
                print("\nAbgebrochen. Tippe 'exit' zum Beenden.")
            except EOFError:
                # Ctrl+D / Dateiende
                print("\nAuf Wiedersehen!")
                break


# =============================================================================
# HAUPTPROGRAMM (direkter Aufruf als Skript)
# =============================================================================

if __name__ == '__main__':
    repl = MathREPL()
    repl.run()
