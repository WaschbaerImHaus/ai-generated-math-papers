# repl.py — Dokumentation

## Übersicht

Das Modul `repl.py` stellt einen interaktiven Read-Eval-Print-Loop (REPL) für
alle mathematischen Module des Projekts bereit. Über einfache Textbefehle können
alle Algorithmen ohne Python-Kenntnisse genutzt werden.

## Verwendung

### Als Skript starten

```bash
python3 src/repl.py
```

### Programmatisch einbinden

```python
from repl import MathREPL

repl = MathREPL()
result = repl.evaluate("prime 97")  # "is_prime(97) = True"
repl.run()  # Interaktive Schleife starten
```

---

## Klasse `MathREPL`

### Konstruktor

```python
repl = MathREPL()
```

Initialisiert alle verfügbaren mathematischen Module:
- `algebra` — Zahlentheorie, Polynome, Gleichungen
- `analysis` — Ableitung, Integration, Grenzwerte
- `linear_algebra` — Vektoren, Matrizen, Eigenwerte
- `fourier` — FFT, DFT
- `proof_theory` — Goldbach, Collatz, Primzahlen
- `statistics_math` — Statistik

---

### `evaluate(command: str) -> str`

Wertet einen Befehlsstring aus und gibt das Ergebnis zurück.

**Rückgabewert:** Formatierter Ergebnis-String.

**Sonderfall:** `evaluate("exit")` gibt `"__EXIT__"` zurück.

---

## Verfügbare Befehle

| Befehl | Syntax | Beispiel | Ergebnis |
|--------|--------|---------|---------|
| `prime` | `prime N` | `prime 97` | `is_prime(97) = True` |
| `factor` | `factor N` | `factor 360` | `prime_factorization(360) = {2: 3, 3: 2, 5: 1}` |
| `gcd` | `gcd A B` | `gcd 48 18` | `gcd(48, 18) = 6` |
| `lcm` | `lcm A B` | `lcm 12 18` | `lcm(12, 18) = 36` |
| `solve` | `solve A B C` | `solve 1 -5 6` | Lösungen von $x^2 - 5x + 6 = 0$ |
| `derive` | `derive FUNC at X` | `derive sin at 0` | $\sin'(0) \approx 1.0$ |
| `integrate` | `integrate F A B` | `integrate sin 0 3.14` | $\int_0^\pi \sin(x)\,dx \approx 2$ |
| `limit` | `limit EXPR VAR PT` | `limit sin(x)/x x 0` | $\lim_{x\to 0}\frac{\sin x}{x} = 1$ |
| `eigenvalues` | `eigenvalues JSON` | `eigenvalues [[1,2],[3,4]]` | Eigenwerte der Matrix |
| `fft` | `fft LIST` | `fft [1,2,3,4,5,6,7,8]` | FFT-Amplituden |
| `goldbach` | `goldbach N` | `goldbach 100` | Goldbach-Zerlegung |
| `collatz` | `collatz N` | `collatz 27` | Stoppzeit = 111 |
| `phi` | `phi N` | `phi 12` | $\varphi(12) = 4$ |
| `modinv` | `modinv A M` | `modinv 3 7` | modulares Inverses |
| `sqrt` | `sqrt N` | `sqrt 2` | $\sqrt{2} \approx 1.41421356$ |
| `help` | `help [BEFEHL]` | `help gcd` | Hilfetext |
| `exit`/`quit` | `exit` | `exit` | REPL beenden |

---

### Unterstützte Funktionen für `derive` und `integrate`

| Name | Funktion |
|------|----------|
| `sin` | $\sin(x)$ |
| `cos` | $\cos(x)$ |
| `tan` | $\tan(x)$ |
| `exp` | $e^x$ |
| `log` | $\ln(x)$ |
| `sqrt` | $\sqrt{x}$ |
| `sq` | $x^2$ |
| `id` | $x$ |
| `abs` | $|x|$ |

---

## Hilfsfunktionen

### `is_jupyter() -> bool`

Erkennt ob das Programm in einem Jupyter-Notebook läuft.

```python
from repl import is_jupyter
if is_jupyter():
    print("Jupyter-Modus aktiv")
```

**Implementierung:** Prüft ob `get_ipython().__class__.__name__ == 'ZMQInteractiveShell'`.

---

### `create_jupyter_notebook(output_path: str) -> None`

Erstellt eine `.ipynb`-Datei mit Demo-Zellen für alle Module.

```python
from repl import create_jupyter_notebook
create_jupyter_notebook("/tmp/demo.ipynb")
```

Das Notebook enthält Zellen für:
- Setup und Imports
- Algebra-Demo (Primzahlen, GGT, quadratische Gleichungen)
- Analysis-Demo (Ableitung, Integration, Grenzwert)
- Lineare Algebra-Demo (Vektoren, Eigenwerte)
- Fourier-Transformation
- LaTeX-Export

---

## Interaktive Session — Beispiel

```
$ python3 src/repl.py
============================================================
  specialist-maths REPL v1.0
  Tippe 'help' für Befehle, 'exit' zum Beenden
============================================================
math> prime 97
is_prime(97) = True

math> gcd 48 18
gcd(48, 18) = 6

math> solve 1 -5 6
solve_quadratic(1.0x² + -5.0x + 6.0 = 0) = [3.0, 2.0]

math> eigenvalues [[4,1],[2,3]]
Eigenwerte = ['5.000000', '2.000000']

math> fft [1,2,3,4]
FFT([1, 2, 3, 4]) =
  Amplituden: ['10.0000', '2.8284', '2.0000', '2.8284']

math> goldbach 100
goldbach(100) = 3 + 97

math> exit
Auf Wiedersehen!
```

---

## Fehlerbehandlung

| Situation | Verhalten |
|-----------|-----------|
| Unbekannter Befehl | Fehlermeldung mit Hinweis auf `help` |
| Falsches Argument | Fehlermeldung mit Syntax-Hinweis |
| Leere Eingabe | Leerer String `""` |
| Python-Exception | Fehlermeldung mit Exception-Typ |
| `exit` / `quit` | Rückgabe `"__EXIT__"` |

---

*Letzte Änderung: 2026-03-09 | Autor: Kurt Ingwer*
