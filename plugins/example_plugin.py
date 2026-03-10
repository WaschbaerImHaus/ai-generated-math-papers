"""
@file example_plugin.py
@brief Beispiel-Plugin: Elementare Zahlentheorie-Erweiterungen.
@description
    Dieses Plugin demonstriert die Plugin-Schnittstelle des specialist-maths
    Systems. Es stellt einige einfache zahlentheoretische Funktionen bereit,
    die über das Plugin-System dynamisch geladen werden können.

    ## Enthaltene Funktionen

    - ``perfect_number_check``: Prüft ob eine Zahl eine vollkommene Zahl ist
    - ``digit_sum``: Berechnet die Quersumme einer natürlichen Zahl
    - ``collatz_length``: Länge der Collatz-Folge für eine gegebene Zahl
    - ``abundant_numbers``: Gibt alle überschüssigen Zahlen bis n zurück

@author Kurt Ingwer
@lastModified 2026-03-10
@version 1.0.0
"""

# Plugin-Name: wird von discover_plugins() zur Registrierung verwendet
PLUGIN_NAME = "example_plugin"

# Plugin-Metadaten
PLUGIN_VERSION = "1.0.0"
PLUGIN_DESCRIPTION = "Elementare zahlentheoretische Hilfsfunktionen"
PLUGIN_AUTHOR = "Kurt Ingwer"


def setup() -> None:
    """
    @brief Initialisierungsfunktion des Beispiel-Plugins.
    @description Keine besonderen Initialisierungen notwendig.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    pass


def perfect_number_check(n: int) -> bool:
    """
    @brief Prüft, ob n eine vollkommene Zahl ist.
    @description
        Eine vollkommene Zahl ist gleich der Summe ihrer echten Teiler.

        Beispiele: 6 = 1+2+3, 28 = 1+2+4+7+14

        Mathematisch: n ist vollkommen, falls σ(n) = 2n,
        wobei σ die Teilersummenfunktion ist.

    @param n  Zu prüfende natürliche Zahl (n ≥ 1).
    @return   True wenn n vollkommen ist, sonst False.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    if n < 2:
        return False
    # Summe aller echten Teiler (Teiler außer n selbst) berechnen
    divisor_sum = sum(d for d in range(1, n) if n % d == 0)
    return divisor_sum == n


def digit_sum(n: int) -> int:
    """
    @brief Berechnet die Quersumme einer natürlichen Zahl.
    @description
        Die Quersumme ist die Summe aller Ziffern der Dezimaldarstellung.

        Beispiel: digit_sum(1234) = 1+2+3+4 = 10

    @param n  Natürliche Zahl (nicht-negativ).
    @return   Quersumme der Zahl.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Betrag verwenden, um negative Zahlen zu unterstützen
    return sum(int(d) for d in str(abs(n)))


def collatz_length(n: int) -> int:
    """
    @brief Berechnet die Länge der Collatz-Folge für n.
    @description
        Die Collatz-Folge startet bei n und wendet wiederholt an:
        - Falls n gerade: n → n/2
        - Falls n ungerade: n → 3n + 1

        Die Folge endet, wenn n = 1 erreicht ist. Die Länge ist die
        Anzahl der Schritte bis einschließlich 1.

    @param n  Startwert (n ≥ 1).
    @return   Anzahl der Schritte bis n = 1.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    if n < 1:
        raise ValueError(f"Collatz-Folge erfordert n ≥ 1, erhalten: {n}")

    steps = 0
    current = n
    # Folge bis zum Zielwert 1 iterieren
    while current != 1:
        if current % 2 == 0:
            # Gerader Fall: halbieren
            current //= 2
        else:
            # Ungerader Fall: 3n + 1
            current = 3 * current + 1
        steps += 1

    return steps


def abundant_numbers(limit: int) -> list:
    """
    @brief Gibt alle überschüssigen Zahlen bis limit zurück.
    @description
        Eine überschüssige (abundante) Zahl ist eine Zahl, bei der die Summe
        ihrer echten Teiler größer als die Zahl selbst ist.

        Beispiel: 12 ist abundant, da 1+2+3+4+6 = 16 > 12.

    @param limit  Obere Grenze (inklusiv).
    @return       Liste aller abundanten Zahlen im Bereich [1, limit].
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    result = []
    for n in range(1, limit + 1):
        # Summe echter Teiler berechnen
        s = sum(d for d in range(1, n) if n % d == 0)
        if s > n:
            result.append(n)
    return result
