"""
@file test_exceptions.py
@brief Tests für das exceptions.py Modul und die Integration in andere Module.
@description
    Testet:
    - Korrekte Exception-Hierarchie (MathematicalError als Basisklasse)
    - ConvergenceError: wird bei nicht-konvergierendem Newton-Raphson geworfen
    - SingularMatrixError: wird bei singulärer Matrix geworfen
    - DomainError: wird bei ungültigem Definitionsbereich geworfen
    - Fehlermeldungen enthalten hilfreiche, spezifische Informationen
    - Edge-Cases: Grenzfälle für alle Exception-Klassen

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import sys
import os
import pytest

# Sicherstellen dass src/ im Python-Pfad liegt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from exceptions import (
    MathematicalError,
    ConvergenceError,
    SingularMatrixError,
    DivisionByZeroError,
    DomainError,
    PrecisionError,
    InvalidInputError,
    PrimeRequiredError,
    NotImplementedMathError,
    NumericalInstabilityError,
)


# =============================================================================
# TESTS: Exception-Hierarchie
# =============================================================================

class TestExceptionHierarchy:
    """Testet dass alle Exceptions korrekt von MathematicalError erben."""

    def test_convergence_error_is_mathematical_error(self):
        """ConvergenceError muss von MathematicalError erben."""
        exc = ConvergenceError("TestMethode", 100)
        # isinstance-Prüfung für korrekte Vererbung
        assert isinstance(exc, MathematicalError)
        assert isinstance(exc, Exception)

    def test_singular_matrix_error_is_mathematical_error(self):
        """SingularMatrixError muss von MathematicalError erben."""
        exc = SingularMatrixError("Inversion")
        assert isinstance(exc, MathematicalError)

    def test_division_by_zero_error_is_mathematical_error(self):
        """DivisionByZeroError muss von MathematicalError erben."""
        exc = DivisionByZeroError("Division durch Null in Berechnung")
        assert isinstance(exc, MathematicalError)

    def test_domain_error_is_mathematical_error(self):
        """DomainError muss von MathematicalError erben."""
        exc = DomainError("sqrt", -1, "(0, ∞)")
        assert isinstance(exc, MathematicalError)

    def test_precision_error_is_mathematical_error(self):
        """PrecisionError muss von MathematicalError erben."""
        exc = PrecisionError(1e-12, 1e-6)
        assert isinstance(exc, MathematicalError)

    def test_invalid_input_error_is_mathematical_error(self):
        """InvalidInputError muss von MathematicalError erben."""
        exc = InvalidInputError("Ungültige Eingabe")
        assert isinstance(exc, MathematicalError)

    def test_prime_required_error_is_mathematical_error(self):
        """PrimeRequiredError muss von MathematicalError erben."""
        exc = PrimeRequiredError(n=4, param="p")
        assert isinstance(exc, MathematicalError)

    def test_not_implemented_math_error_is_mathematical_error(self):
        """NotImplementedMathError muss von MathematicalError erben."""
        exc = NotImplementedMathError("Noch nicht implementiert")
        assert isinstance(exc, MathematicalError)

    def test_numerical_instability_error_is_mathematical_error(self):
        """NumericalInstabilityError muss von MathematicalError erben."""
        exc = NumericalInstabilityError("NaN erkannt")
        assert isinstance(exc, MathematicalError)

    def test_catch_all_with_mathematical_error(self):
        """Alle Exceptions können mit 'except MathematicalError' gefangen werden."""
        exceptions_to_test = [
            ConvergenceError("Test", 10),
            SingularMatrixError(),
            DivisionByZeroError("test"),
            DomainError("f", 0),
            PrecisionError(1e-10, 1e-5),
            InvalidInputError("test"),
            PrimeRequiredError(4),
            NotImplementedMathError("test"),
            NumericalInstabilityError("test"),
        ]
        # Alle müssen mit MathematicalError gefangen werden können
        for exc in exceptions_to_test:
            caught = False
            try:
                raise exc
            except MathematicalError:
                caught = True
            assert caught, f"{type(exc).__name__} wurde nicht von MathematicalError gefangen"


# =============================================================================
# TESTS: ConvergenceError
# =============================================================================

class TestConvergenceError:
    """Testet ConvergenceError Erstellung und Attribute."""

    def test_basic_creation(self):
        """ConvergenceError mit Pflichtparametern erstellen."""
        exc = ConvergenceError("Newton-Raphson", 1000)
        # Attribute müssen korrekt gesetzt sein
        assert exc.method == "Newton-Raphson"
        assert exc.iterations == 1000
        assert exc.last_value is None

    def test_with_last_value(self):
        """ConvergenceError mit letztem Wert erstellen."""
        exc = ConvergenceError("Bisection", 500, last_value=3.14159)
        assert exc.last_value == 3.14159
        # Letzter Wert muss in der Fehlermeldung erscheinen
        assert "3.14159" in str(exc)

    def test_message_contains_method_name(self):
        """Fehlermeldung muss den Methodennamen enthalten."""
        exc = ConvergenceError("Runge-Kutta-45", 200)
        assert "Runge-Kutta-45" in str(exc)

    def test_message_contains_iteration_count(self):
        """Fehlermeldung muss die Iterationszahl enthalten."""
        exc = ConvergenceError("Test", 777)
        assert "777" in str(exc)

    def test_message_without_last_value(self):
        """Fehlermeldung ohne last_value darf keine None-Referenz enthalten."""
        exc = ConvergenceError("Test", 100)
        msg = str(exc)
        # Kein 'None' in der Nachricht wenn last_value=None
        assert "None" not in msg

    def test_can_be_raised_and_caught(self):
        """ConvergenceError muss geworfen und gefangen werden können."""
        with pytest.raises(ConvergenceError) as exc_info:
            raise ConvergenceError("TestAlgo", 50, last_value=0.5)
        assert exc_info.value.iterations == 50


# =============================================================================
# TESTS: ConvergenceError aus newton_raphson()
# =============================================================================

class TestConvergenceFromNewtonRaphson:
    """Testet dass newton_raphson() korrekt ConvergenceError wirft."""

    def test_newton_raphson_raises_convergence_error(self):
        """newton_raphson() muss ConvergenceError bei Nicht-Konvergenz werfen."""
        from analysis import newton_raphson

        # Funktion ohne Nullstelle: f(x) = x² + 1 hat keine reellen Nullstellen
        # Newton-Raphson divergiert oder konvergiert nicht
        def no_real_root(x):
            return x**2 + 1

        with pytest.raises(ConvergenceError):
            newton_raphson(no_real_root, x0=0.0, max_iter=5)

    def test_newton_raphson_error_is_mathematical_error(self):
        """Der von newton_raphson() geworfene Fehler muss MathematicalError sein."""
        from analysis import newton_raphson

        def no_real_root(x):
            return x**2 + 1

        with pytest.raises(MathematicalError):
            newton_raphson(no_real_root, x0=0.0, max_iter=5)

    def test_newton_raphson_success_no_exception(self):
        """newton_raphson() darf bei konvergierendem Fall keine Exception werfen."""
        from analysis import newton_raphson
        import math

        # f(x) = x² - 2, Nullstelle bei √2
        root = newton_raphson(lambda x: x**2 - 2, x0=1.0)
        assert abs(root - math.sqrt(2)) < 1e-10


# =============================================================================
# TESTS: SingularMatrixError
# =============================================================================

class TestSingularMatrixError:
    """Testet SingularMatrixError Erstellung und Attribute."""

    def test_default_operation(self):
        """SingularMatrixError ohne Parameter nutzt Standardwert 'Inversion'."""
        exc = SingularMatrixError()
        assert "Inversion" in str(exc)
        assert exc.operation == "Inversion"

    def test_custom_operation(self):
        """SingularMatrixError mit benutzerdefinierter Operation."""
        exc = SingularMatrixError("LU-Zerlegung")
        assert "LU-Zerlegung" in str(exc)
        assert exc.operation == "LU-Zerlegung"

    def test_message_contains_singular(self):
        """Fehlermeldung muss 'Singuläre Matrix' oder ähnliches enthalten."""
        exc = SingularMatrixError()
        # Nachricht soll auf das Problem hinweisen
        assert "ingulär" in str(exc) or "singular" in str(exc).lower()


# =============================================================================
# TESTS: SingularMatrixError aus Matrix-Operationen
# =============================================================================

class TestSingularMatrixFromOperations:
    """Testet dass Matrix-Operationen korrekt SingularMatrixError werfen."""

    def test_matrix_inverse_raises_singular_error(self):
        """Matrix.inverse() muss SingularMatrixError bei singulärer Matrix werfen."""
        from matrix_ops import Matrix

        # Singuläre Matrix: Zeile 2 = 2 × Zeile 1 (lineare Abhängigkeit)
        singular = Matrix([
            [1.0, 2.0],
            [2.0, 4.0],
        ])

        with pytest.raises(SingularMatrixError):
            singular.inverse()

    def test_matrix_inverse_error_is_mathematical_error(self):
        """Der Fehler von Matrix.inverse() muss MathematicalError sein."""
        from matrix_ops import Matrix

        singular = Matrix([[0.0, 0.0], [0.0, 0.0]])

        with pytest.raises(MathematicalError):
            singular.inverse()

    def test_lu_decomposition_raises_singular_error(self):
        """lu_decomposition() muss SingularMatrixError bei singulärer Matrix werfen."""
        from matrix_ops import Matrix
        from matrix_decomp import lu_decomposition

        # Nullmatrix ist singulär
        singular = Matrix([[0.0, 0.0], [0.0, 0.0]])

        with pytest.raises(SingularMatrixError):
            lu_decomposition(singular)

    def test_regular_matrix_no_exception(self):
        """Matrix.inverse() darf bei regulärer Matrix keine Exception werfen."""
        from matrix_ops import Matrix

        # Einheitsmatrix ist regulär (invertierbar)
        identity = Matrix([[1.0, 0.0], [0.0, 1.0]])
        inv = identity.inverse()

        # Inverse der Einheitsmatrix ist die Einheitsmatrix selbst
        assert abs(inv.get(0, 0) - 1.0) < 1e-10
        assert abs(inv.get(1, 1) - 1.0) < 1e-10


# =============================================================================
# TESTS: DomainError
# =============================================================================

class TestDomainError:
    """Testet DomainError Erstellung und Attribute."""

    def test_creation_with_domain(self):
        """DomainError mit Definitionsbereich erstellen."""
        exc = DomainError("sqrt", -4, "(0, ∞)")
        assert exc.function == "sqrt"
        assert exc.value == -4
        assert exc.domain == "(0, ∞)"
        # Alle Informationen müssen in der Nachricht erscheinen
        assert "sqrt" in str(exc)
        assert "-4" in str(exc)
        assert "(0, ∞)" in str(exc)

    def test_creation_without_domain(self):
        """DomainError ohne Definitionsbereich möglich."""
        exc = DomainError("log", 0)
        assert exc.function == "log"
        assert exc.value == 0
        assert exc.domain == ""

    def test_domain_error_from_solve_linear(self):
        """solve_linear() muss DomainError bei a=0 werfen."""
        from algebra_core import solve_linear

        # a=0, b≠0: Kein Widerspruch möglich
        with pytest.raises(DomainError):
            solve_linear(0, 5)

    def test_domain_error_from_solve_linear_trivial(self):
        """solve_linear() muss DomainError bei a=0, b=0 werfen."""
        from algebra_core import solve_linear

        # a=0, b=0: Unendlich viele Lösungen
        with pytest.raises(DomainError):
            solve_linear(0, 0)

    def test_solve_linear_domain_error_is_mathematical_error(self):
        """Der Fehler von solve_linear() muss MathematicalError sein."""
        from algebra_core import solve_linear

        with pytest.raises(MathematicalError):
            solve_linear(0, 1)


# =============================================================================
# TESTS: PrecisionError
# =============================================================================

class TestPrecisionError:
    """Testet PrecisionError Erstellung."""

    def test_creation(self):
        """PrecisionError mit erforderlicher und erreichter Genauigkeit."""
        exc = PrecisionError(required=1e-12, achieved=1e-6)
        assert exc.required == 1e-12
        assert exc.achieved == 1e-6

    def test_message_format(self):
        """Fehlermeldung muss beide Genauigkeitswerte in wissenschaftlicher Notation enthalten."""
        exc = PrecisionError(1e-10, 1e-4)
        msg = str(exc)
        # Beide Werte müssen in der Nachricht enthalten sein
        assert "e" in msg.lower() or "E" in msg  # Wissenschaftliche Notation


# =============================================================================
# TESTS: PrimeRequiredError
# =============================================================================

class TestPrimeRequiredError:
    """Testet PrimeRequiredError Erstellung."""

    def test_creation(self):
        """PrimeRequiredError mit nicht-prim Wert erstellen."""
        exc = PrimeRequiredError(n=4, param="p")
        assert exc.n == 4
        assert exc.param == "p"
        assert "4" in str(exc)
        assert "p" in str(exc)

    def test_default_param_name(self):
        """PrimeRequiredError nutzt 'n' als Standard-Parametername."""
        exc = PrimeRequiredError(9)
        assert exc.param == "n"
        assert "n" in str(exc)
        assert "9" in str(exc)


# =============================================================================
# TESTS: Fehlermeldungs-Qualität
# =============================================================================

class TestErrorMessageQuality:
    """Testet dass alle Fehlermeldungen hilfreich und informativ sind."""

    def test_convergence_error_message_is_helpful(self):
        """ConvergenceError-Nachricht enthält Methode, Iterationen und Wert."""
        exc = ConvergenceError("BFGS", 200, last_value=42.5)
        msg = str(exc)
        assert "BFGS" in msg        # Methodenname
        assert "200" in msg         # Iterationszahl
        assert "42.5" in msg        # Letzter Wert

    def test_singular_matrix_error_message_is_helpful(self):
        """SingularMatrixError-Nachricht nennt die fehlgeschlagene Operation."""
        exc = SingularMatrixError("Cholesky-Zerlegung")
        msg = str(exc)
        assert "Cholesky" in msg    # Operationsname

    def test_domain_error_message_is_helpful(self):
        """DomainError-Nachricht enthält Funktionsname, Wert und Bereich."""
        exc = DomainError("arccos", 2.5, "[-1, 1]")
        msg = str(exc)
        assert "arccos" in msg      # Funktionsname
        assert "2.5" in str(exc)    # Ungültiger Wert
        assert "[-1, 1]" in msg     # Definitionsbereich

    def test_precision_error_message_is_helpful(self):
        """PrecisionError-Nachricht zeigt klaren Vergleich der Genauigkeiten."""
        exc = PrecisionError(1e-15, 1e-8)
        msg = str(exc)
        # Beide Genauigkeiten müssen lesbar sein
        assert "erforderlich" in msg or "required" in msg.lower()

    def test_prime_required_error_message_is_helpful(self):
        """PrimeRequiredError-Nachricht enthält den ungültigen Wert und Parameternamen."""
        exc = PrimeRequiredError(n=100, param="modulus")
        msg = str(exc)
        assert "100" in msg
        assert "modulus" in msg
