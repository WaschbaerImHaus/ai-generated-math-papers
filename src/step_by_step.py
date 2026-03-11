"""
@file step_by_step.py
@brief Schritt-für-Schritt-Erklärungen mathematischer Algorithmen.
@description
    Gibt detaillierte, lehrreiche Erklärungen zu Berechnungsschritten aus.
    Ideal für Lernzwecke: zeigt nicht nur das Ergebnis, sondern den gesamten
    Lösungsweg mit mathematischen Erklärungen.

    Jede Funktion gibt ein dict zurück mit:
        'steps'     → Liste der Einzelschritte als dict
        'result'    → Endergebnis
        'method'    → Name des Algorithmus
        'converged' → (bei iterativen Verfahren) bool

    Verwendung:
        from step_by_step import newton_raphson_steps, format_steps_text
        result = newton_raphson_steps(lambda x: x**2 - 2, 1.5)
        print(format_steps_text(result))

@author Michael Fuhrmann
@date 2026-03-10
@lastModified 2026-03-10
"""

import sys
import os
import math
import copy

# Elternverzeichnis (src/) im Suchpfad aufnehmen
sys.path.insert(0, os.path.dirname(__file__))


# ===========================================================================
# HILFSFUNKTIONEN (intern)
# ===========================================================================

def _numerical_derivative(f, x: float, h: float = 1e-7) -> float:
    """
    @brief Zentrale Differenzenformel für numerische Ableitung.
    @description
        f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
        Fehler O(h²), deutlich besser als Vorwärtsdifferenz O(h).

    @param f: Zu differenzierende Funktion
    @param x: Stelle, an der abgeleitet wird
    @param h: Schrittweite (Standard: 1e-7 für guten Kompromiss Rundung/Trunkierung)
    @return Numerische Ableitung f'(x)
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    return (f(x + h) - f(x - h)) / (2 * h)


def _matrix_to_str(matrix: list) -> str:
    """
    @brief Wandelt eine 2D-Liste in einen formatierten Matrix-String um.
    @description
        Jede Zeile wird als [ a  b  c ] dargestellt, Zahlen auf 4 Dezimalstellen.

    @param matrix: 2D-Liste
    @return Mehrzeiliger String der Matrix
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    rows = []
    for row in matrix:
        vals = []
        for v in row:
            if isinstance(v, float):
                vals.append(f"{v:9.4f}")
            else:
                vals.append(f"{v:9}")
        rows.append("[ " + "  ".join(vals) + " ]")
    return "\n".join(rows)


def _augmented_to_str(matrix: list, b: list) -> str:
    """
    @brief Stellt die erweiterte Koeffizientenmatrix [A|b] als String dar.
    @description
        Verwendet in der Gauss-Elimination zur Visualisierung.

    @param matrix: Koeffizientenmatrix A
    @param b: Rechte-Seite-Vektor b
    @return Formatierten String von [A|b]
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    rows = []
    for i, row in enumerate(matrix):
        vals = [f"{v:9.4f}" if isinstance(v, float) else f"{v:9}" for v in row]
        b_val = f"{b[i]:9.4f}" if isinstance(b[i], float) else f"{b[i]:9}"
        rows.append("[ " + "  ".join(vals) + " | " + b_val + " ]")
    return "\n".join(rows)


# ===========================================================================
# SCHRITT-FÜR-SCHRITT-FUNKTIONEN
# ===========================================================================

def newton_raphson_steps(
    f,
    x0: float,
    tol: float = 1e-10,
    max_iter: int = 50
) -> dict:
    """
    @brief Newton-Raphson-Verfahren mit vollständiger Schritt-Dokumentation.
    @description
        Das Newton-Raphson-Verfahren findet Nullstellen einer Funktion f iterativ:
            x_{n+1} = x_n - f(x_n) / f'(x_n)

        Geometrische Interpretation:
        - In jedem Schritt wird die Tangente an f im Punkt (x_n, f(x_n)) gezogen
        - Die Nullstelle dieser Tangente wird als nächste Näherung verwendet
        - Konvergenz: quadratisch nahe der Nullstelle (sehr schnell!)

        KaTeX-Formel: x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}

    @param f: Funktion, deren Nullstelle gesucht wird (callable)
    @param x0: Startwert (Anfangsnäherung)
    @param tol: Konvergenztoleranz – Abbruch wenn |f(x)| < tol
    @param max_iter: Maximale Anzahl Iterationen
    @return Dict mit 'steps', 'result', 'converged', 'iterations', 'method'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    steps = []
    x = float(x0)
    converged = False

    for i in range(1, max_iter + 1):
        fx = f(x)
        dfx = _numerical_derivative(f, x)

        # Division durch Null abfangen (Tangente horizontal)
        if abs(dfx) < 1e-15:
            steps.append({
                'iteration': i,
                'x': x,
                'fx': fx,
                'dfx': dfx,
                'x_new': None,
                'explanation': (
                    f"Schritt {i}: f'(x) ≈ 0 – Die Tangente ist waagerecht. "
                    f"Newton-Raphson kann nicht fortgesetzt werden."
                ),
                'error': True
            })
            break

        # Newton-Schritt berechnen: x_{n+1} = x_n - f(x_n)/f'(x_n)
        x_new = x - fx / dfx

        # Menschenlesbare Erklärung des Schritts aufbauen
        explanation = (
            f"Schritt {i}: "
            f"x{i} = x{i-1} - f(x{i-1})/f'(x{i-1}) "
            f"= {x:.6f} - ({fx:.6f})/({dfx:.6f}) "
            f"= {x_new:.6f}"
        )

        steps.append({
            'iteration': i,
            'x': x,
            'fx': fx,
            'dfx': dfx,
            'x_new': x_new,
            'residual': abs(fx),
            'explanation': explanation
        })

        # Konvergenzkritierium: Funktion nahe genug an Null?
        if abs(fx) < tol:
            converged = True
            break

        x = x_new

    return {
        'steps': steps,
        'result': x,
        'converged': converged,
        'iterations': len(steps),
        'method': 'Newton-Raphson',
        'formula': r'x_{n+1} = x_n - \frac{f(x_n)}{f\'(x_n)}'
    }


def bisection_steps(
    f,
    a: float,
    b: float,
    tol: float = 1e-10,
    max_iter: int = 50
) -> dict:
    """
    @brief Bisektionsverfahren mit vollständiger Schritt-Dokumentation.
    @description
        Das Bisektionsverfahren (auch: Intervallhalbierungsverfahren) findet
        eine Nullstelle von f im Intervall [a, b], wobei f(a)·f(b) < 0.

        Jeder Schritt halbiert das Intervall:
        1. Berechne Mittelpunkt m = (a+b)/2
        2. Prüfe Vorzeichen von f(m)
        3. Wähle die Hälfte, in der die Nullstelle liegt

        Konvergenz: linear, Fehler halbiert sich pro Schritt.
        KaTeX-Formel: m = \frac{a+b}{2}

    @param f: Funktion (callable)
    @param a: Linke Intervallgrenze (f(a) und f(b) müssen verschiedene Vorzeichen haben)
    @param b: Rechte Intervallgrenze
    @param tol: Toleranz – Abbruch wenn |b-a| < tol
    @param max_iter: Maximale Anzahl Iterationen
    @return Dict mit 'steps', 'result', 'converged', 'iterations', 'method'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    steps = []

    # Vorzeichenprüfung: Bisektion erfordert Vorzeichenwechsel
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        return {
            'steps': [],
            'result': None,
            'converged': False,
            'iterations': 0,
            'method': 'Bisektion',
            'error': f"f(a)={fa:.4f} und f(b)={fb:.4f} haben gleiches Vorzeichen – keine Nullstelle garantiert."
        }

    converged = False
    midpoint = a

    for i in range(1, max_iter + 1):
        # Intervallmitte berechnen
        midpoint = (a + b) / 2.0
        fm = f(midpoint)
        interval_width = b - a

        # Erklärung des aktuellen Schritts
        if fm * fa < 0:
            new_a, new_b = a, midpoint
            side = "linke Hälfte [a, m]"
        elif fm * fb < 0:
            new_a, new_b = midpoint, b
            side = "rechte Hälfte [m, b]"
        else:
            # Exakte Nullstelle gefunden
            new_a, new_b = midpoint, midpoint
            side = "exakte Nullstelle"

        explanation = (
            f"Schritt {i}: m = ({a:.6f} + {b:.6f}) / 2 = {midpoint:.6f}, "
            f"f(m) = {fm:.6f}, "
            f"Vorzeichenwechsel in {side} → neues Intervall [{new_a:.6f}, {new_b:.6f}]"
        )

        steps.append({
            'iteration': i,
            'a': a,
            'b': b,
            'midpoint': midpoint,
            'f_midpoint': fm,
            'interval_width': interval_width,
            'new_a': new_a,
            'new_b': new_b,
            'explanation': explanation
        })

        # Konvergenz erreicht wenn Intervall klein genug oder Nullstelle exakt
        if abs(fm) < tol or interval_width / 2 < tol:
            converged = True
            break

        # Intervall aktualisieren
        a, b = new_a, new_b
        fa = f(a)
        fb = f(b)

    return {
        'steps': steps,
        'result': midpoint,
        'converged': converged,
        'iterations': len(steps),
        'method': 'Bisektion (Intervallhalbierung)',
        'formula': r'm = \frac{a + b}{2}'
    }


def gauss_elimination_steps(matrix: list[list[float]], b: list[float]) -> dict:
    """
    @brief Gauss-Elimination mit Schritt-Dokumentation (Vorwärts + Rückwärts).
    @description
        Löst das lineare Gleichungssystem Ax = b durch:
        1. Vorwärts-Elimination: Bringt A auf obere Dreiecksform
           (mit Teilpivotsuche für numerische Stabilität)
        2. Rückwärts-Substitution: Löst das Dreieckssystem

        Teilpivotsuche: In jeder Spalte wird die Zeile mit dem betragsmäßig
        größten Element als Pivotzeile gewählt → vermeidet Division durch sehr
        kleine Zahlen (numerische Instabilität).

        KaTeX-Formel: A \\mathbf{x} = \\mathbf{b}

    @param matrix: Koeffizientenmatrix A als 2D-Liste (n×n)
    @param b: Rechte-Seite-Vektor als Liste (Länge n)
    @return Dict mit 'steps', 'result', 'method'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    steps = []
    n = len(matrix)

    # Tiefe Kopie, damit Originaldaten nicht verändert werden
    A = [list(map(float, row)) for row in matrix]
    bv = list(map(float, b))

    # -----------------------------------------------------------------------
    # PHASE 1: Vorwärts-Elimination (auf obere Dreiecksform bringen)
    # -----------------------------------------------------------------------
    for col in range(n):
        # Pivot-Suche: Zeile mit größtem Betrag in Spalte col (ab Zeile col)
        pivot_row = col
        max_val = abs(A[col][col])
        for row in range(col + 1, n):
            if abs(A[row][col]) > max_val:
                max_val = abs(A[row][col])
                pivot_row = row

        # Zeilen tauschen falls nötig
        if pivot_row != col:
            A[col], A[pivot_row] = A[pivot_row], A[col]
            bv[col], bv[pivot_row] = bv[pivot_row], bv[col]

            steps.append({
                'phase': 'Vorwärts-Elimination',
                'pivot_col': col,
                'pivot_row': col,
                'operation': f"Zeile {col+1} ↔ Zeile {pivot_row+1} (Pivotsuche)",
                'matrix': _augmented_to_str(A, bv),
                'explanation': (
                    f"Teilpivotsuche in Spalte {col+1}: "
                    f"Zeile {pivot_row+1} hat größten Betrag |{max_val:.4f}|. "
                    f"Tausch mit Zeile {col+1} für numerische Stabilität."
                )
            })

        # Pivot-Element
        pivot = A[col][col]
        if abs(pivot) < 1e-14:
            # Singuläre Matrix – Gleichungssystem nicht eindeutig lösbar
            steps.append({
                'phase': 'Vorwärts-Elimination',
                'pivot_col': col,
                'operation': f"FEHLER: Pivot ≈ 0 in Spalte {col+1}",
                'matrix': _augmented_to_str(A, bv),
                'explanation': "Matrix ist singulär – das Gleichungssystem hat keine eindeutige Lösung."
            })
            return {
                'steps': steps,
                'result': None,
                'method': 'Gauss-Elimination',
                'error': 'Singuläre Matrix'
            }

        # Alle Zeilen unterhalb der Pivotzeile eliminieren
        for row in range(col + 1, n):
            # Multiplikator für diese Zeile berechnen
            factor = A[row][col] / pivot

            if abs(factor) < 1e-15:
                # Multiplikator ist Null → Zeile bereits in Ordnung
                continue

            # Zeilenoperation: Zeile[row] -= factor * Zeile[col]
            for k in range(col, n):
                A[row][k] -= factor * A[col][k]
            bv[row] -= factor * bv[col]

            steps.append({
                'phase': 'Vorwärts-Elimination',
                'pivot_col': col,
                'pivot_row': col,
                'target_row': row,
                'factor': factor,
                'operation': f"Z{row+1} = Z{row+1} - ({factor:.4f}) × Z{col+1}",
                'matrix': _augmented_to_str(A, bv),
                'explanation': (
                    f"Eliminiere Element in Zeile {row+1}, Spalte {col+1}: "
                    f"Multiplikator m = {A[row][col] + factor * A[col][col]:.4f} / {pivot:.4f} = {factor:.4f}. "
                    f"Z{row+1} ← Z{row+1} - {factor:.4f} · Z{col+1}"
                )
            })

    # -----------------------------------------------------------------------
    # PHASE 2: Rückwärts-Substitution
    # -----------------------------------------------------------------------
    x = [0.0] * n

    steps.append({
        'phase': 'Rückwärts-Substitution (Start)',
        'operation': 'Obere Dreiecksmatrix lösen',
        'matrix': _augmented_to_str(A, bv),
        'explanation': (
            "Die Matrix ist jetzt in oberer Dreiecksform. "
            "Rückwärts-Substitution: x_n aus letzter Zeile, dann aufwärts."
        )
    })

    for row in range(n - 1, -1, -1):
        # Rechte Seite der bereits bekannten Lösungen abziehen
        rhs = bv[row]
        for k in range(row + 1, n):
            rhs -= A[row][k] * x[k]

        # Aktuelles x[row] berechnen
        x[row] = rhs / A[row][row]

        steps.append({
            'phase': 'Rückwärts-Substitution',
            'variable': row,
            'value': x[row],
            'operation': f"x{row+1} = ({bv[row]:.4f} - Σ) / {A[row][row]:.4f} = {x[row]:.6f}",
            'matrix': _augmented_to_str(A, bv),
            'explanation': (
                f"x{row+1} = {x[row]:.6f} "
                f"(Zeile {row+1}: {A[row][row]:.4f} · x{row+1} = {bv[row]:.4f} - Summe bekannter Terme)"
            )
        })

    return {
        'steps': steps,
        'result': x,
        'method': 'Gauss-Elimination mit Teilpivotsuche',
        'formula': r'A\mathbf{x} = \mathbf{b}'
    }


def euclidean_algorithm_steps(a: int, b: int) -> dict:
    """
    @brief Euklidischer Algorithmus mit vollständiger Schritt-Dokumentation.
    @description
        Der euklidische Algorithmus berechnet den größten gemeinsamen Teiler (ggT)
        von zwei ganzen Zahlen a und b.

        Grundprinzip:
            ggT(a, b) = ggT(b, a mod b)
            Terminiert wenn b = 0, dann ist a der ggT.

        Jeder Schritt hat die Form: a = b·q + r
        wobei q = a // b (ganzzahlige Division) und r = a mod b (Rest).

        KaTeX-Formel: \\gcd(a, b) = \\gcd(b, a \\bmod b)

    @param a: Erste ganze Zahl (positiv)
    @param b: Zweite ganze Zahl (positiv)
    @return Dict mit 'steps', 'result' (ggT), 'method'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    steps = []
    original_a, original_b = a, b

    # Sicherstellen, dass beide Zahlen positiv sind
    a, b = abs(a), abs(b)

    step_num = 1
    while b != 0:
        # Division mit Rest: a = b*q + r
        q = a // b
        r = a % b

        explanation = (
            f"Schritt {step_num}: {a} = {b} × {q} + {r} "
            f"[Division: {a} ÷ {b} = {q} Rest {r}]"
        )

        steps.append({
            'step': step_num,
            'a': a,
            'b': b,
            'quotient': q,
            'remainder': r,
            'operation': f"{a} = {b} × {q} + {r}",
            'explanation': explanation
        })

        # Nächste Iteration: ggT(a,b) = ggT(b, r)
        a, b = b, r
        step_num += 1

    # Letzter Schritt: b = 0, a ist der ggT
    steps.append({
        'step': step_num,
        'a': a,
        'b': 0,
        'quotient': None,
        'remainder': 0,
        'operation': f"b = 0 → ggT = {a}",
        'explanation': (
            f"Schritt {step_num}: b = 0. "
            f"Der ggT von {original_a} und {original_b} ist {a}."
        )
    })

    return {
        'steps': steps,
        'result': a,  # ggT
        'method': 'Euklidischer Algorithmus',
        'formula': r'\gcd(a, b) = \gcd(b, a \bmod b)',
        'input_a': original_a,
        'input_b': original_b
    }


def rsa_steps(p: int, q: int, message: int) -> dict:
    """
    @brief RSA-Kryptographie Schritt-für-Schritt (Schlüsselgenerierung + Ver-/Entschlüsselung).
    @description
        RSA (Rivest-Shamir-Adleman) ist ein asymmetrisches Kryptosystem.
        Sicherheit basiert auf der Schwierigkeit, große Zahlen zu faktorisieren.

        Schritte:
        1. n = p·q (RSA-Modul)
        2. φ(n) = (p-1)·(q-1) (Euler'sche Phi-Funktion)
        3. e wählen: 1 < e < φ(n), ggT(e, φ(n)) = 1
        4. d berechnen: d·e ≡ 1 (mod φ(n)) [modulares Inverses]
        5. Verschlüsseln: c = m^e mod n
        6. Entschlüsseln: m = c^d mod n

        KaTeX-Formeln:
            n = p \\cdot q
            \\phi(n) = (p-1)(q-1)
            c \\equiv m^e \\pmod{n}
            m \\equiv c^d \\pmod{n}

    @param p: Erste Primzahl
    @param q: Zweite Primzahl
    @param message: Zu verschlüsselnde Nachricht (ganzzahlig, < n)
    @return Dict mit allen RSA-Schritten und 'result'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    steps = []

    # -----------------------------------------------------------------------
    # Schritt 1: RSA-Modul berechnen
    # -----------------------------------------------------------------------
    n = p * q
    steps.append({
        'step': 1,
        'title': 'RSA-Modul berechnen',
        'operation': f"n = p × q = {p} × {q} = {n}",
        'explanation': (
            f"Der öffentliche Modul n = {n}. "
            f"Dieser ist das Produkt zweier Primzahlen. "
            f"Nachrichten müssen kleiner als n = {n} sein."
        )
    })

    # -----------------------------------------------------------------------
    # Schritt 2: Euler'sche Phi-Funktion
    # -----------------------------------------------------------------------
    phi_n = (p - 1) * (q - 1)
    steps.append({
        'step': 2,
        'title': "Euler'sche Phi-Funktion",
        'operation': f"φ(n) = (p-1)(q-1) = ({p}-1)({q}-1) = {p-1} × {q-1} = {phi_n}",
        'explanation': (
            f"φ({n}) = {phi_n} gibt an, wie viele Zahlen von 1 bis {n} "
            f"teilerfremd zu {n} sind. "
            f"Für n = p·q gilt: φ(n) = (p-1)(q-1)."
        )
    })

    # -----------------------------------------------------------------------
    # Schritt 3: Öffentlichen Exponenten e wählen
    # -----------------------------------------------------------------------
    # Standard-Kandidaten für e (bevorzugt: 65537, dann kleinere)
    e_candidates = [65537, 257, 17, 5, 3]
    e = None
    for candidate in e_candidates:
        if candidate < phi_n and math.gcd(candidate, phi_n) == 1:
            e = candidate
            break

    # Fallback: Kleinstes e > 1 mit ggT(e, phi_n) = 1
    if e is None:
        for candidate in range(2, phi_n):
            if math.gcd(candidate, phi_n) == 1:
                e = candidate
                break

    steps.append({
        'step': 3,
        'title': 'Öffentlichen Exponenten e wählen',
        'operation': f"e = {e} (gewählt: ggT({e}, {phi_n}) = {math.gcd(e, phi_n)} = 1)",
        'explanation': (
            f"Der öffentliche Exponent e = {e} muss: "
            f"(a) 1 < e < φ(n) = {phi_n} und "
            f"(b) ggT(e, φ(n)) = 1 erfüllen. "
            f"e = {e} erfüllt beides: ggT({e}, {phi_n}) = {math.gcd(e, phi_n)}."
        )
    })

    # -----------------------------------------------------------------------
    # Schritt 4: Privaten Schlüssel d berechnen (mod. Inverses via erw. Euklid)
    # -----------------------------------------------------------------------
    def extended_gcd(a: int, b: int):
        """Erweiterter euklidischer Algorithmus: gibt (gcd, x, y) zurück mit a·x + b·y = gcd."""
        if b == 0:
            return a, 1, 0
        g, x1, y1 = extended_gcd(b, a % b)
        return g, y1, x1 - (a // b) * y1

    _, d, _ = extended_gcd(e, phi_n)
    # d muss positiv sein
    d = d % phi_n

    steps.append({
        'step': 4,
        'title': 'Privaten Schlüssel d berechnen',
        'operation': f"d = e⁻¹ mod φ(n) = {e}⁻¹ mod {phi_n} = {d}",
        'explanation': (
            f"Modulares Inverses von e={e} modulo φ(n)={phi_n}: "
            f"d·e ≡ 1 (mod {phi_n}). "
            f"Berechnung via erweitertem euklidischen Algorithmus: d = {d}. "
            f"Probe: {d} × {e} mod {phi_n} = {(d * e) % phi_n} (= 1 ✓)"
        )
    })

    # -----------------------------------------------------------------------
    # Schritt 5: Schlüsselpaare zusammenfassen
    # -----------------------------------------------------------------------
    steps.append({
        'step': 5,
        'title': 'Schlüsselpaare',
        'operation': f"Öffentlich: (e={e}, n={n})  |  Privat: (d={d}, n={n})",
        'explanation': (
            f"Öffentlicher Schlüssel: (e={e}, n={n}) – kann frei geteilt werden. "
            f"Privater Schlüssel: (d={d}, n={n}) – geheim halten!"
        )
    })

    # -----------------------------------------------------------------------
    # Schritt 6: Verschlüsseln
    # -----------------------------------------------------------------------
    if message >= n:
        message_adjusted = message % n
        steps.append({
            'step': 6,
            'title': 'Nachricht anpassen',
            'operation': f"m = {message} mod {n} = {message_adjusted}",
            'explanation': f"Nachricht muss kleiner als n={n} sein. Angepasst: {message_adjusted}."
        })
        message = message_adjusted

    ciphertext = pow(message, e, n)
    steps.append({
        'step': 6 if message < n else 7,
        'title': 'Verschlüsselung',
        'operation': f"c = m^e mod n = {message}^{e} mod {n} = {ciphertext}",
        'explanation': (
            f"Verschlüsselung: c = m^e mod n = {message}^{e} mod {n} = {ciphertext}. "
            f"Der Geheimtext ist {ciphertext}."
        )
    })

    # -----------------------------------------------------------------------
    # Schritt 7: Entschlüsseln
    # -----------------------------------------------------------------------
    decrypted = pow(ciphertext, d, n)
    steps.append({
        'step': 7,
        'title': 'Entschlüsselung',
        'operation': f"m = c^d mod n = {ciphertext}^{d} mod {n} = {decrypted}",
        'explanation': (
            f"Entschlüsselung: m = c^d mod n = {ciphertext}^{d} mod {n} = {decrypted}. "
            f"{'Erfolgreich! ✓' if decrypted == message else 'Fehler! ✗'} "
            f"Entschlüsselt: {decrypted} (Original: {message})"
        )
    })

    return {
        'steps': steps,
        'result': {
            'n': n,
            'phi_n': phi_n,
            'e': e,
            'd': d,
            'public_key': (e, n),
            'private_key': (d, n),
            'message': message,
            'ciphertext': ciphertext,
            'decrypted': decrypted,
            'success': decrypted == message
        },
        'method': 'RSA-Kryptographie',
        'formula': r'c \\equiv m^e \\pmod{n},\quad m \\equiv c^d \\pmod{n}'
    }


def prime_factorization_steps(n: int) -> dict:
    """
    @brief Primfaktorzerlegung mit vollständiger Schritt-Dokumentation.
    @description
        Zerlegt eine positive ganze Zahl n in ihre Primfaktoren.
        Methode: Trial Division (Probedivision):
        1. Teile durch 2 so oft wie möglich
        2. Teile durch ungerade Zahlen ab 3 (bis sqrt(n))
        3. Falls n > 1 verbleibt: n selbst ist ein Primfaktor

        Beispiel: 360 = 2³ × 3² × 5

        KaTeX-Formel: n = \\prod_{i} p_i^{e_i}

    @param n: Positive ganze Zahl, die zerlegt werden soll
    @return Dict mit 'steps', 'result' (Liste der Primfaktoren), 'method'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    steps = []
    original_n = n
    factors = []
    remaining = n
    step_num = 1

    # -----------------------------------------------------------------------
    # Phase 1: Division durch 2 (einzige gerade Primzahl)
    # -----------------------------------------------------------------------
    if remaining % 2 == 0:
        count = 0
        while remaining % 2 == 0:
            remaining //= 2
            count += 1
            factors.append(2)

            steps.append({
                'step': step_num,
                'divisor': 2,
                'quotient': remaining,
                'operation': f"{remaining * 2} ÷ 2 = {remaining}",
                'explanation': (
                    f"Schritt {step_num}: {remaining * 2} ist durch 2 teilbar. "
                    f"{remaining * 2} ÷ 2 = {remaining}. "
                    f"2 ist ein Primfaktor."
                )
            })
            step_num += 1

    # -----------------------------------------------------------------------
    # Phase 2: Ungerade Faktoren ab 3 testen (bis sqrt(remaining))
    # -----------------------------------------------------------------------
    divisor = 3
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            remaining //= divisor
            factors.append(divisor)

            steps.append({
                'step': step_num,
                'divisor': divisor,
                'quotient': remaining,
                'operation': f"{remaining * divisor} ÷ {divisor} = {remaining}",
                'explanation': (
                    f"Schritt {step_num}: {remaining * divisor} ist durch {divisor} teilbar. "
                    f"{remaining * divisor} ÷ {divisor} = {remaining}. "
                    f"{divisor} ist ein Primfaktor."
                )
            })
            step_num += 1

        divisor += 2  # Nur ungerade Zahlen testen

    # -----------------------------------------------------------------------
    # Phase 3: Falls remaining > 1, ist es selbst prim
    # -----------------------------------------------------------------------
    if remaining > 1:
        factors.append(remaining)
        steps.append({
            'step': step_num,
            'divisor': remaining,
            'quotient': 1,
            'operation': f"{remaining} ist prim",
            'explanation': (
                f"Schritt {step_num}: Verbleibender Wert {remaining} > 1 "
                f"und hatte keinen Teiler bis √{original_n} ≈ {math.isqrt(original_n)}. "
                f"Also ist {remaining} selbst ein Primfaktor."
            )
        })

    # Darstellung als Produkt mit Exponenten
    factor_counts = {}
    for f in factors:
        factor_counts[f] = factor_counts.get(f, 0) + 1

    product_str = " × ".join(
        f"{p}^{e}" if e > 1 else str(p)
        for p, e in sorted(factor_counts.items())
    )

    # Abschluss-Schritt
    steps.append({
        'step': step_num + (1 if remaining > 1 else 0),
        'divisor': None,
        'quotient': 1,
        'operation': f"{original_n} = {product_str}",
        'explanation': (
            f"Ergebnis: {original_n} = {product_str}. "
            f"Alle Primfaktoren: {sorted(set(factors))}"
        )
    })

    return {
        'steps': steps,
        'result': factors,
        'factors_unique': sorted(set(factors)),
        'factor_counts': factor_counts,
        'product_str': product_str,
        'method': 'Primfaktorzerlegung (Trial Division)',
        'formula': r'n = \prod_{i} p_i^{e_i}'
    }


def lu_decomposition_steps(matrix: list[list[float]]) -> dict:
    """
    @brief LU-Zerlegung (Doolittle) mit vollständiger Schritt-Dokumentation.
    @description
        Zerlegt eine quadratische Matrix A in das Produkt:
            PA = LU
        wobei:
        - P: Permutationsmatrix (Zeilentausch, Teilpivotsuche)
        - L: Untere Dreiecksmatrix (Lower) mit 1en auf der Diagonale
        - U: Obere Dreiecksmatrix (Upper)

        Doolittle-Algorithmus:
        - U[i,j] = A[i,j] - Σ L[i,k]·U[k,j]  (für k < i)
        - L[i,j] = (A[i,j] - Σ L[i,k]·U[k,j]) / U[j,j]  (für i > j)

        KaTeX-Formel: PA = LU

    @param matrix: Quadratische Matrix als 2D-Liste
    @return Dict mit 'steps', 'result' (L, U, P), 'method'
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    steps = []
    n = len(matrix)

    # Tiefe Kopie für Berechnungen
    A = [list(map(float, row)) for row in matrix]

    # Identitätsmatrizen für L und Permutationsvektor
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [row[:] for row in A]
    perm = list(range(n))  # Permutationsvektor

    steps.append({
        'phase': 'Initialisierung',
        'operation': 'Ausgangszustand',
        'L_matrix': _matrix_to_str(L),
        'U_matrix': _matrix_to_str(U),
        'explanation': (
            f"Starte LU-Zerlegung einer {n}×{n}-Matrix. "
            f"L = Einheitsmatrix, U = Kopie von A. "
            f"Ziel: A = L·U (Doolittle-Algorithmus)."
        )
    })

    # -----------------------------------------------------------------------
    # Doolittle-Algorithmus mit Teilpivotsuche
    # -----------------------------------------------------------------------
    for k in range(n):
        # Pivot-Suche in Spalte k
        pivot_row = k
        max_val = abs(U[k][k])
        for i in range(k + 1, n):
            if abs(U[i][k]) > max_val:
                max_val = abs(U[i][k])
                pivot_row = i

        # Zeilen tauschen falls nötig
        if pivot_row != k:
            U[k], U[pivot_row] = U[pivot_row], U[k]
            perm[k], perm[pivot_row] = perm[pivot_row], perm[k]

            # In L nur den bereits berechneten Teil tauschen
            for j in range(k):
                L[k][j], L[pivot_row][j] = L[pivot_row][j], L[k][j]

            steps.append({
                'phase': f'Spalte {k+1}: Pivotsuche',
                'operation': f"Zeile {k+1} ↔ Zeile {pivot_row+1}",
                'L_matrix': _matrix_to_str(L),
                'U_matrix': _matrix_to_str(U),
                'explanation': (
                    f"Teilpivotsuche: Zeile {pivot_row+1} hat größten Betrag "
                    f"|{max_val:.4f}| in Spalte {k+1}. Tausch mit Zeile {k+1}."
                )
            })

        # Multiplikatoren berechnen und Elimination durchführen
        for i in range(k + 1, n):
            if abs(U[k][k]) < 1e-14:
                # Singuläre Matrix
                break

            # Multiplikator (wird in L gespeichert)
            factor = U[i][k] / U[k][k]
            L[i][k] = factor

            # Zeile i in U aktualisieren
            for j in range(k, n):
                U[i][j] -= factor * U[k][j]

            steps.append({
                'phase': f'Spalte {k+1}: Elimination',
                'operation': f"L[{i+1},{k+1}] = {factor:.4f}, U: Z{i+1} -= {factor:.4f} × Z{k+1}",
                'L_matrix': _matrix_to_str(L),
                'U_matrix': _matrix_to_str(U),
                'explanation': (
                    f"Multiplikator für Zeile {i+1}: "
                    f"l_{{i,k}} = u_{{i,k}} / u_{{k,k}} = {U[i][k] + factor * U[k][k]:.4f} / {U[k][k]:.4f} = {factor:.4f}. "
                    f"In L gespeichert. U-Zeile {i+1} aktualisiert."
                )
            })

    # Permutationsmatrix aus Permutationsvektor erstellen
    P = [[1.0 if perm[i] == j else 0.0 for j in range(n)] for i in range(n)]

    steps.append({
        'phase': 'Ergebnis',
        'operation': 'LU-Zerlegung abgeschlossen: PA = LU',
        'L_matrix': _matrix_to_str(L),
        'U_matrix': _matrix_to_str(U),
        'explanation': (
            "Die Zerlegung PA = LU ist abgeschlossen. "
            "L hat 1en auf der Diagonale (Einheits-Dreiecksmatrix). "
            "U ist die obere Dreiecksmatrix."
        )
    })

    return {
        'steps': steps,
        'result': {'L': L, 'U': U, 'P': P, 'perm': perm},
        'method': 'LU-Zerlegung (Doolittle mit Teilpivotsuche)',
        'formula': r'PA = LU'
    }


# ===========================================================================
# FORMATIERUNGSFUNKTIONEN
# ===========================================================================

def format_steps_text(steps_dict: dict) -> str:
    """
    @brief Formatiert ein Schritt-für-Schritt-Ergebnis als lesbaren Konsolen-Text.
    @description
        Wandelt das von den step-Funktionen zurückgegebene Dict in einen
        menschenlesbaren mehrzeiligen String um.

        Enthält:
        - Methodenname als Überschrift
        - Jeden Schritt nummeriert mit Operation und Erklärung
        - Endergebnis am Schluss

    @param steps_dict: Ausgabe-Dict einer step-Funktion
    @return Formatierter String für Konsole
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    lines = []
    method = steps_dict.get('method', 'Algorithmus')
    formula = steps_dict.get('formula', '')

    # Überschrift
    lines.append("=" * 60)
    lines.append(f"  {method}")
    if formula:
        lines.append(f"  Formel: {formula}")
    lines.append("=" * 60)
    lines.append("")

    # Alle Schritte ausgeben
    for step in steps_dict.get('steps', []):
        # Schritt-Header (Iteration oder Step-Nummer)
        num = step.get('iteration') or step.get('step') or step.get('phase', '')
        title = step.get('title', '')

        if title:
            lines.append(f"--- Schritt {num}: {title} ---")
        elif num:
            lines.append(f"--- Schritt {num} ---")

        # Operation (kompakte Darstellung)
        operation = step.get('operation', '')
        if operation:
            lines.append(f"  Operation: {operation}")

        # Erklärung (ausführlich)
        explanation = step.get('explanation', '')
        if explanation:
            lines.append(f"  Erklärung: {explanation}")

        # Matrix-Darstellung falls vorhanden
        for key in ('matrix', 'L_matrix', 'U_matrix'):
            if key in step:
                label = {'matrix': 'Matrix', 'L_matrix': 'L', 'U_matrix': 'U'}.get(key, key)
                lines.append(f"  {label}:")
                for mat_line in step[key].split('\n'):
                    lines.append(f"    {mat_line}")

        lines.append("")

    # Endergebnis
    result = steps_dict.get('result')
    lines.append("-" * 60)
    if isinstance(result, dict):
        lines.append("  Ergebnis:")
        for k, v in result.items():
            if isinstance(v, float):
                lines.append(f"    {k} = {v:.6f}")
            else:
                lines.append(f"    {k} = {v}")
    elif isinstance(result, list):
        lines.append(f"  Ergebnis: {result}")
    elif isinstance(result, float):
        lines.append(f"  Ergebnis: {result:.8f}")
    else:
        lines.append(f"  Ergebnis: {result}")

    # Konvergenz-Info falls vorhanden
    if 'converged' in steps_dict:
        status = "konvergiert" if steps_dict['converged'] else "NICHT konvergiert"
        lines.append(f"  Status: {status} nach {steps_dict.get('iterations', '?')} Schritten")

    lines.append("=" * 60)
    return "\n".join(lines)


def format_steps_html(steps_dict: dict) -> str:
    """
    @brief Formatiert ein Schritt-für-Schritt-Ergebnis als HTML für Web-Darstellung.
    @description
        Erstellt semantisches HTML mit:
        - <div class="step"> für jeden Schritt
        - <code> für Operationen und Formeln
        - KaTeX-kompatible Formel-Delimitier ($...$)
        - Klassen für CSS-Styling und JavaScript-Animationen

    @param steps_dict: Ausgabe-Dict einer step-Funktion
    @return HTML-String mit vollständigem Schritt-Container
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    method = steps_dict.get('method', 'Algorithmus')
    formula = steps_dict.get('formula', '')
    steps = steps_dict.get('steps', [])

    html_parts = []
    html_parts.append(f'<div class="steps-container">')
    html_parts.append(f'  <div class="steps-header">')
    html_parts.append(f'    <h3 class="steps-title">{method}</h3>')

    # KaTeX-Formel anzeigen
    if formula:
        html_parts.append(f'    <div class="steps-formula">$${ formula}$$</div>')
    html_parts.append(f'  </div>')

    html_parts.append(f'  <div class="steps-list">')

    for idx, step in enumerate(steps):
        num = step.get('iteration') or step.get('step') or (idx + 1)
        title = step.get('title', '')
        operation = step.get('operation', '')
        explanation = step.get('explanation', '')
        phase = step.get('phase', '')
        is_error = step.get('error', False)

        # CSS-Klassen für Animationen und Fehler-Markierung
        step_class = "step step-error" if is_error else "step"

        html_parts.append(f'    <div class="{step_class}" data-step="{idx}">')

        # Schritt-Nummer und Titel
        step_label = title or phase or f"Schritt {num}"
        html_parts.append(
            f'      <div class="step-header">'
            f'<span class="step-number">{num}</span>'
            f'<span class="step-title">{step_label}</span>'
            f'</div>'
        )

        # Operation oder Erklärung als Code (Fallback: explanation)
        code_content = operation or explanation
        if code_content:
            # Sonderzeichen für HTML escapen
            op_escaped = code_content.replace('<', '&lt;').replace('>', '&gt;')
            html_parts.append(
                f'      <div class="step-operation">'
                f'<code>{op_escaped}</code>'
                f'</div>'
            )

        # Erklärung als Text
        if explanation:
            exp_escaped = explanation.replace('<', '&lt;').replace('>', '&gt;')
            html_parts.append(
                f'      <div class="step-explanation">{exp_escaped}</div>'
            )

        # Matrix-Darstellung falls vorhanden
        for key in ('matrix', 'L_matrix', 'U_matrix'):
            if key in step:
                label = {'matrix': 'Matrix', 'L_matrix': 'L', 'U_matrix': 'U'}.get(key, key)
                mat_escaped = step[key].replace('<', '&lt;').replace('>', '&gt;')
                html_parts.append(
                    f'      <div class="step-matrix">'
                    f'<span class="matrix-label">{label}:</span>'
                    f'<pre><code>{mat_escaped}</code></pre>'
                    f'</div>'
                )

        html_parts.append(f'    </div>')  # Ende .step

    html_parts.append(f'  </div>')  # Ende .steps-list

    # Ergebnis-Box
    result = steps_dict.get('result')
    converged = steps_dict.get('converged')
    iterations = steps_dict.get('iterations', '?')

    html_parts.append(f'  <div class="steps-result">')
    html_parts.append(f'    <h4>Ergebnis</h4>')

    if isinstance(result, dict):
        html_parts.append(f'    <dl class="result-dict">')
        for k, v in result.items():
            v_str = f"{v:.6f}" if isinstance(v, float) else str(v)
            html_parts.append(f'      <dt>{k}</dt><dd><code>{v_str}</code></dd>')
        html_parts.append(f'    </dl>')
    elif isinstance(result, list):
        result_escaped = str(result).replace('<', '&lt;').replace('>', '&gt;')
        html_parts.append(f'    <code class="result-value">{result_escaped}</code>')
    elif isinstance(result, float):
        html_parts.append(f'    <code class="result-value">{result:.8f}</code>')
    elif result is not None:
        html_parts.append(f'    <code class="result-value">{result}</code>')

    # Konvergenz-Status
    if converged is not None:
        status_class = "converged" if converged else "not-converged"
        status_text = f"Konvergiert nach {iterations} Schritten" if converged else f"Nicht konvergiert ({iterations} Schritte)"
        html_parts.append(
            f'    <div class="convergence-status {status_class}">{status_text}</div>'
        )

    html_parts.append(f'  </div>')  # Ende .steps-result
    html_parts.append(f'</div>')  # Ende .steps-container

    return "\n".join(html_parts)
