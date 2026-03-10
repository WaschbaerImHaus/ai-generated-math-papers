"""
@file exceptions.py
@brief Spezifische Ausnahmen für mathematische Fehler.
@description
    Definiert eine Hierarchie von Exception-Klassen für alle
    mathematischen Fehlerfälle im specialist-maths Projekt.
    Ermöglicht präzise Fehlerbehandlung statt generischer ValueError/RuntimeError.

    Hierarchie:
        MathematicalError (Basisklasse)
        ├── ConvergenceError        – Algorithmus konvergiert nicht
        ├── SingularMatrixError     – Matrix ist singulär (nicht invertierbar)
        ├── DivisionByZeroError     – Division durch Null
        ├── DomainError             – Eingabe außerhalb Definitionsbereich
        ├── PrecisionError          – Numerische Genauigkeit unterschritten
        ├── InvalidInputError       – Ungültige Eingabe
        ├── PrimeRequiredError      – Primzahl als Argument erforderlich
        ├── NotImplementedMathError – Math. Fall noch nicht implementiert
        └── NumericalInstabilityError – Numerische Instabilität (NaN, Inf)

@author Kurt Ingwer
@lastModified 2026-03-10
"""


class MathematicalError(Exception):
    """
    @brief Basisklasse aller mathematischen Ausnahmen im specialist-maths Projekt.
    @description
        Alle anderen mathematischen Fehlerklassen erben von dieser Klasse.
        Ermöglicht das Abfangen aller math. Fehler mit einem einzigen except-Block:

            try:
                ...
            except MathematicalError as e:
                print(f"Mathematischer Fehler: {e}")

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    pass


class ConvergenceError(MathematicalError):
    """
    @brief Ausnahme wenn ein iterativer Algorithmus nicht konvergiert.
    @description
        Wird geworfen von: newton_raphson(), bisection(), runge_kutta45() usw.
        wenn die maximale Iterationszahl ohne Konvergenz erreicht wird.

    @example
        raise ConvergenceError("Newton-Raphson", 1000, last_value=3.14)
        # Ergibt: "Newton-Raphson konvergiert nicht nach 1000 Iterationen, letzter Wert: 3.14"

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, method: str, iterations: int, last_value: float = None):
        """
        @brief Erstellt eine ConvergenceError-Ausnahme.
        @param method Name des Algorithmus/Verfahrens (z.B. "Newton-Raphson").
        @param iterations Anzahl der durchgeführten Iterationen bis zum Abbruch.
        @param last_value Letzter berechneter Wert (optional, für Debugging).
        @date 2026-03-10
        """
        # Attribute für programmatischen Zugriff speichern
        self.method = method
        self.iterations = iterations
        self.last_value = last_value

        # Fehlermeldung aufbauen (mit optionalem letzten Wert)
        msg = f"{method} konvergiert nicht nach {iterations} Iterationen"
        if last_value is not None:
            msg += f", letzter Wert: {last_value}"
        super().__init__(msg)


class SingularMatrixError(MathematicalError):
    """
    @brief Ausnahme wenn eine singuläre (nicht invertierbare) Matrix vorliegt.
    @description
        Eine Matrix ist singulär wenn ihre Determinante 0 ist.
        Wird geworfen bei: Matrix.inverse(), lu_decomposition(), Matrix.solve() usw.

        Singuläre Matrizen haben keinen vollen Rang und können für LGS
        entweder keine oder unendlich viele Lösungen haben.

    @example
        raise SingularMatrixError("LU-Zerlegung")
        # Ergibt: "Singuläre Matrix: LU-Zerlegung nicht möglich"

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, operation: str = "Inversion"):
        """
        @brief Erstellt eine SingularMatrixError-Ausnahme.
        @param operation Beschreibung der fehlgeschlagenen Operation.
        @date 2026-03-10
        """
        self.operation = operation
        super().__init__(f"Singuläre Matrix: {operation} nicht möglich")


class DivisionByZeroError(MathematicalError):
    """
    @brief Ausnahme bei mathematischer Division durch Null.
    @description
        Unterscheidet sich von Pythons eingebautem ZeroDivisionError dadurch,
        dass sie explizit mathematische Kontexte anzeigt (z.B. Modulo-Rechnung,
        Nenner in rationalen Ausdrücken).

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    pass


class DomainError(MathematicalError):
    """
    @brief Ausnahme wenn ein Eingabewert außerhalb des Definitionsbereichs liegt.
    @description
        Wird geworfen wenn z.B.:
        - sqrt(-1) real berechnet werden soll
        - log(0) oder log(-x) aufgerufen wird
        - Eine Funktion nur für positive Zahlen definiert ist und ein negativer Wert übergeben wird

    @example
        raise DomainError("log", -1, domain="(0, ∞)")
        # Ergibt: "log: Wert -1 liegt außerhalb des Definitionsbereichs (0, ∞)"

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, function: str, value, domain: str = ""):
        """
        @brief Erstellt eine DomainError-Ausnahme.
        @param function Name der betroffenen Funktion.
        @param value Der ungültige Eingabewert.
        @param domain Beschreibung des Definitionsbereichs (optional).
        @date 2026-03-10
        """
        self.function = function
        self.value = value
        self.domain = domain

        # Fehlermeldung mit optionalem Definitionsbereich
        msg = f"{function}: Wert {value} liegt außerhalb des Definitionsbereichs"
        if domain:
            msg += f" {domain}"
        super().__init__(msg)


class PrecisionError(MathematicalError):
    """
    @brief Ausnahme wenn numerische Genauigkeit die Anforderung unterschreitet.
    @description
        Wird geworfen wenn ein Algorithmus zwar terminiert, aber die
        erzielte Genauigkeit nicht ausreicht (z.B. bei schlechter Kondition).

    @example
        raise PrecisionError(required=1e-12, achieved=1e-6)
        # Ergibt: "Genauigkeit: erforderlich 1.00e-12, erreicht 1.00e-06"

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, required: float, achieved: float):
        """
        @brief Erstellt eine PrecisionError-Ausnahme.
        @param required Erforderliche Genauigkeit (Toleranz).
        @param achieved Tatsächlich erreichte Genauigkeit.
        @date 2026-03-10
        """
        self.required = required
        self.achieved = achieved
        super().__init__(
            f"Genauigkeit: erforderlich {required:.2e}, erreicht {achieved:.2e}"
        )


class InvalidInputError(MathematicalError):
    """
    @brief Ausnahme bei ungültigen Eingaben für mathematische Funktionen.
    @description
        Allgemeine Ausnahme für ungültige Parameter, die nicht in die
        spezifischeren Kategorien (DomainError, PrimeRequiredError) fallen.

        Beispiele:
        - Negative Anzahl von Termen
        - Leere Listen wo Elemente erwartet werden
        - Typfehler in mathematischen Funktionen

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    pass


class PrimeRequiredError(MathematicalError):
    """
    @brief Ausnahme wenn eine Primzahl als Argument erwartet wird.
    @description
        Wird geworfen von Funktionen die explizit Primzahlen erwarten,
        z.B. p_adic_valuation(), rsa_keygen() mit kleinen Primzahlen usw.

    @example
        raise PrimeRequiredError(n=4, param="p")
        # Ergibt: "Parameter p=4 muss eine Primzahl sein"

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """

    def __init__(self, n: int, param: str = "n"):
        """
        @brief Erstellt eine PrimeRequiredError-Ausnahme.
        @param n Der übergebene (nicht-prime) Wert.
        @param param Name des Parameters (z.B. "p", "q", "n").
        @date 2026-03-10
        """
        self.n = n
        self.param = param
        super().__init__(f"Parameter {param}={n} muss eine Primzahl sein")


class NotImplementedMathError(MathematicalError):
    """
    @brief Ausnahme wenn ein mathematischer Fall noch nicht implementiert ist.
    @description
        Ähnlich wie Pythons NotImplementedError, aber speziell für
        mathematische Randfälle (z.B. bestimmte Sonderfälle einer Funktion,
        die theoretisch lösbar aber noch nicht implementiert sind).

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    pass


class NumericalInstabilityError(MathematicalError):
    """
    @brief Ausnahme bei erkannter numerischer Instabilität.
    @description
        Wird geworfen wenn ein Algorithmus NaN, Inf oder extreme Werte
        produziert, die auf numerische Instabilität hinweisen.

        Ursachen:
        - Überlauf (Overflow) bei großen Zahlen
        - Auslöschung bei fast-gleichen Zahlen
        - Schlecht konditionierte Matrizen oder Funktionen

    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    pass
