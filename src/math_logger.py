"""
@file math_logger.py
@brief Logging-System für mathematische Berechnungen.
@description
    Protokolliert Berechnungsschritte für Debugging und Lernzwecke.
    Unterstützt verschiedene Log-Level und Ausgabeformate.

    Verwendung:
        from math_logger import MathLogger
        logger = MathLogger(level="DEBUG")
        logger.step("Newton-Raphson", iteration=1, x=2.5, fx=0.25)
        logger.result("Nullstelle", value=3.0)

    Log-Level (aufsteigend):
        DEBUG   → alle Berechnungsschritte (verbose)
        INFO    → Ergebnisse und Fortschritt
        WARNING → Konvergenz-Warnungen, ungewöhnliche Zustände
        ERROR   → Fehler, die die Berechnung abbrechen

@author Michael Fuhrmann
@date 2026-03-10
@lastModified 2026-03-10
"""

import logging
import time
from datetime import datetime
from typing import Any
from pathlib import Path


# ===========================================================================
# HAUPTKLASSE: MathLogger
# ===========================================================================

class MathLogger:
    """
    @brief Logger für mathematische Berechnungen.
    @description
        Kapselt Python-logging für die spezifischen Anforderungen des
        Mathematik-Spezialisten. Bietet formatierte Ausgaben für:
        - Berechnungsschritte (Iteration, Werte)
        - Endergebnisse
        - Konvergenz-Warnungen
        - Matrix-Umformungsschritte
        - Ausführungszeiten

        Ausgabe: Konsole und/oder Datei unter logs/

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """

    # Abbildung von String-Level auf Python-logging-Konstanten
    _LEVEL_MAP = {
        "DEBUG":   logging.DEBUG,
        "INFO":    logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR":   logging.ERROR,
    }

    def __init__(
        self,
        name: str = "specialist-maths",
        level: str = "INFO",
        log_to_file: bool = False,
        log_dir: str = None
    ):
        """
        @brief Initialisiert den Logger.
        @description
            Erstellt einen Python-Logger mit konfigurierbarem Level und optionalem
            Datei-Handler. Der Logger ist sofort einsatzbereit nach dem Erzeugen.

        @param name: Name des Loggers (erscheint im Log-Prefix)
        @param level: Log-Level als String: "DEBUG"|"INFO"|"WARNING"|"ERROR"
        @param log_to_file: Wenn True, werden Logs auch in Datei geschrieben
        @param log_dir: Verzeichnis für Log-Dateien (Standard: projekt/logs/)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self.name = name

        # Python-Logger mit eindeutigem Namen erzeugen
        # Eindeutigkeit verhindert Konflikte bei mehreren MathLogger-Instanzen
        self._logger = logging.getLogger(f"mathlogger.{name}.{id(self)}")
        self._logger.setLevel(self._LEVEL_MAP.get(level.upper(), logging.INFO))

        # Verhindern, dass Log-Meldungen an Parent-Logger weitergegeben werden
        # (würde sonst in root-Logger landen und doppelt ausgegeben)
        self._logger.propagate = False

        # Alle existierenden Handler entfernen (z.B. bei Re-Initialisierung)
        self._logger.handlers.clear()

        # Konsolen-Handler immer hinzufügen
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(self._build_formatter())
        self._logger.addHandler(console_handler)

        # Optional: Datei-Handler hinzufügen
        if log_to_file:
            log_path = self._resolve_log_dir(log_dir)
            log_file = log_path / f"{name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_handler.setFormatter(self._build_formatter())
            self._logger.addHandler(file_handler)

    # ---------------------------------------------------------------------------
    # PRIVATE HILFSMETHODEN
    # ---------------------------------------------------------------------------

    def _build_formatter(self) -> logging.Formatter:
        """
        @brief Erzeugt einen einheitlichen Formatter für alle Handler.
        @description
            Format: [LEVEL] HH:MM:SS | Nachricht
            Zeitstempel hilft beim Nachvollziehen der Ausführungsreihenfolge.

        @return logging.Formatter-Instanz
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        return logging.Formatter(
            fmt="[%(levelname)-7s] %(asctime)s | %(message)s",
            datefmt="%H:%M:%S"
        )

    def _resolve_log_dir(self, log_dir: str | None) -> Path:
        """
        @brief Bestimmt und erstellt das Log-Verzeichnis.
        @description
            Standard: Projektverzeichnis/logs/ (zwei Ebenen über webapp/src/).
            Erstellt das Verzeichnis falls es noch nicht existiert.

        @param log_dir: Überschreibt den Standardpfad (oder None für Standard)
        @return Path-Objekt des Log-Verzeichnisses
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        if log_dir:
            path = Path(log_dir)
        else:
            # __file__ ist src/math_logger.py → zwei Ebenen hoch = Projektverzeichnis
            project_root = Path(__file__).parent.parent
            path = project_root / "logs"

        # Verzeichnis erzeugen falls nicht vorhanden (existiert-ok)
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _format_kwargs(self, kwargs: dict) -> str:
        """
        @brief Formatiert Keyword-Argumente als lesbaren String.
        @description
            Zahlen werden mit 4 Dezimalstellen formatiert für einheitliche Darstellung.
            Strings werden unverändert ausgegeben.

        @param kwargs: Dictionary mit Schrittwerten
        @return Formatierter String "key1=val1, key2=val2, ..."
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        parts = []
        for key, val in kwargs.items():
            # Gleitkommazahlen einheitlich formatieren
            if isinstance(val, float):
                parts.append(f"{key}={val:.4f}")
            else:
                parts.append(f"{key}={val}")
        return ", ".join(parts)

    def _format_matrix(self, matrix: list) -> str:
        """
        @brief Formatiert eine Matrix als lesbares Gitter.
        @description
            Jede Zeile wird eingerückt und die Werte werden ausgerichtet.
            Beispiel:
                [ 2.0000  1.0000 ]
                [ 0.0000  3.0000 ]

        @param matrix: 2D-Liste mit Matrixwerten
        @return Formatierter String
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        if not matrix:
            return "[]"

        # Maximale Breite ermitteln für gleichmäßige Ausrichtung
        formatted_rows = []
        for row in matrix:
            formatted_rows.append([f"{val:8.4f}" if isinstance(val, float) else f"{val:8}" for val in row])

        lines = []
        for row in formatted_rows:
            lines.append("  [ " + "  ".join(row) + " ]")

        return "\n" + "\n".join(lines)

    # ---------------------------------------------------------------------------
    # ÖFFENTLICHE LOGGING-METHODEN
    # ---------------------------------------------------------------------------

    def step(self, algorithm: str, **kwargs) -> None:
        """
        @brief Protokolliert einen Berechnungsschritt.
        @description
            Wird bei jeder Iteration eines iterativen Verfahrens aufgerufen.
            Gibt auf DEBUG-Level aus, da Schritte sehr zahlreich sein können.

            Format: [DEBUG] Newton-Raphson | Schritt: iteration=3, x=1.5000, fx=0.2500

        @param algorithm: Name des Algorithmus (z.B. "Newton-Raphson")
        @param kwargs: Schrittparameter (z.B. iteration=1, x=2.5, fx=0.25)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        params_str = self._format_kwargs(kwargs)
        self._logger.debug(f"{algorithm} | Schritt: {params_str}")

    def result(self, description: str, value: Any, unit: str = "") -> None:
        """
        @brief Protokolliert ein Endergebnis einer Berechnung.
        @description
            Wird nach Abschluss einer Berechnung aufgerufen (INFO-Level).
            Format: [INFO ] Ergebnis: description = value unit

        @param description: Bezeichnung des Ergebnisses (z.B. "Nullstelle")
        @param value: Ergebniswert (beliebiger Typ)
        @param unit: Optionale Einheit (z.B. "rad", "m/s")
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        # Einheit mit Leerzeichen voranstellen falls angegeben
        unit_str = f" {unit}" if unit else ""

        # Gleitkommazahlen auf 6 Dezimalstellen formatieren
        if isinstance(value, float):
            value_str = f"{value:.6f}"
        else:
            value_str = str(value)

        self._logger.info(f"Ergebnis: {description} = {value_str}{unit_str}")

    def convergence_warning(self, algorithm: str, iteration: int, residual: float) -> None:
        """
        @brief Warnung bei langsamer Konvergenz eines iterativen Verfahrens.
        @description
            Wird ausgegeben wenn das Residuum nach vielen Iterationen noch groß ist.
            Format: [WARNING] algorithm: Langsame Konvergenz bei Iteration N, Residuum=1.2345e-03

        @param algorithm: Name des iterativen Verfahrens
        @param iteration: Aktuelle Iterationsnummer
        @param residual: Aktuelles Residuum (|f(x)| oder Fehler)
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self._logger.warning(
            f"{algorithm}: Langsame Konvergenz bei Iteration {iteration}, "
            f"Residuum={residual:.4e}"
        )

    def matrix_step(
        self,
        operation: str,
        matrix_before: list,
        matrix_after: list,
        description: str = ""
    ) -> None:
        """
        @brief Protokolliert einen Matrix-Umformungsschritt.
        @description
            Stellt die Matrix vor und nach der Operation lesbar dar.
            Nützlich für Gauss-Elimination, LU-Zerlegung, etc.

        @param operation: Beschreibung der Operation (z.B. "Zeile 2 -= 2 * Zeile 1")
        @param matrix_before: Matrix vor der Operation (2D-Liste)
        @param matrix_after: Matrix nach der Operation (2D-Liste)
        @param description: Optionale Erklärung der Operation
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        desc_str = f" ({description})" if description else ""
        before_str = self._format_matrix(matrix_before)
        after_str = self._format_matrix(matrix_after)

        self._logger.debug(
            f"Matrix-Op: {operation}{desc_str}\n"
            f"  Vorher:{before_str}\n"
            f"  Nachher:{after_str}"
        )

    def timing(self, func_name: str, elapsed: float) -> None:
        """
        @brief Protokolliert die Ausführungszeit einer Funktion.
        @description
            Gibt auf INFO-Level aus, da Timing für Leistungsanalyse relevant ist.
            Format: [INFO ] Timing: func_name dauerte 0.0123s

        @param func_name: Name der gemessenen Funktion
        @param elapsed: Verstrichene Zeit in Sekunden
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self._logger.info(f"Timing: {func_name} dauerte {elapsed:.4f}s")

    def section(self, title: str) -> None:
        """
        @brief Gibt eine Trennlinie mit Titel aus für besser lesbare Logs.
        @description
            Erzeugt optische Trennung zwischen Berechnungsabschnitten.
            Format: [INFO ] ===== TITEL =====

        @param title: Abschnittsbezeichnung
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        separator = "=" * 50
        self._logger.info(f"{separator}")
        self._logger.info(f"  {title.upper()}")
        self._logger.info(f"{separator}")

    def warning(self, message: str) -> None:
        """
        @brief Gibt eine allgemeine Warnung auf WARNING-Level aus.
        @description
            Öffentliche Schnittstelle für Warnmeldungen, die nicht in die
            spezialisierten Warners (convergence_warning, condition_number_warning)
            passen. Wird z.B. von _safe_parse() verwendet, wenn parse_expr()
            scheitert und auf sp.sympify() zurückgegriffen werden muss.

            Format: [WARNING] <message>

        @param message: Warntext als beliebiger String.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Interne Python-logging-Instanz direkt nutzen
        self._logger.warning(message)

    def set_level(self, level: str) -> None:
        """
        @brief Ändert den Log-Level zur Laufzeit.
        @description
            Ermöglicht das dynamische Aktivieren/Deaktivieren von Debug-Ausgaben
            ohne den Logger neu zu erstellen.

        @param level: Neuer Level als String: "DEBUG"|"INFO"|"WARNING"|"ERROR"
        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        new_level = self._LEVEL_MAP.get(level.upper(), logging.INFO)
        self._logger.setLevel(new_level)
        # Alle Handler auf den neuen Level setzen
        for handler in self._logger.handlers:
            handler.setLevel(new_level)

    def condition_number_warning(
        self,
        matrix_name: str,
        cond: float,
        threshold: float = 1e10
    ) -> None:
        """
        @brief Warnung bei schlecht konditionierter Matrix.
        @description
            Gibt eine WARNING-Meldung aus, wenn die Konditionszahl κ(A) den
            Schwellenwert überschreitet. Eine hohe Konditionszahl bedeutet,
            dass kleine Fehler in den Eingabedaten zu großen Fehlern im
            Ergebnis führen können.

            Konditionszahl κ(A) = ||A|| · ||A⁻¹|| = σ_max / σ_min

            Faustregel: Wenn κ ≈ 10^k, verliert man k Dezimalstellen Genauigkeit.
            Bei double precision (16 Dezimalstellen) und κ = 10^10 bleiben nur
            noch ~6 gültige Dezimalstellen übrig.

        @param matrix_name: Name/Bezeichnung der Matrix für den Log.
        @param cond: Berechnete Konditionszahl.
        @param threshold: Schwellenwert (Standard: 1e10).
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Nur warnen wenn Konditionszahl den Schwellenwert überschreitet
        if cond > threshold:
            self._logger.warning(
                f"Schlecht konditionierte Matrix '{matrix_name}': "
                f"κ(A) = {cond:.3e} > {threshold:.1e}. "
                f"Numerische Ergebnisse können ungenau sein."
            )

    def stability_report(self, operations: list) -> str:
        """
        @brief Gibt eine Zusammenfassung aller Konditionszahlen aus.
        @description
            Analysiert eine Liste von Operationen mit zugehörigen Konditionszahlen
            und gibt einen formatierten Bericht aus.

            Format jedes Eintrags: {"name": str, "cond": float, "operation": str}

        @param operations: Liste von Dicts mit Konditionszahlen-Einträgen.
        @return: Formatierter Bericht als String.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Bericht-Zeilen aufbauen
        lines = ["=== Stabilitätsbericht ==="]

        # Alle Operationen nach Konditionszahl aufsteigend sortieren
        sorted_ops = sorted(
            [op for op in operations if isinstance(op, dict) and 'cond' in op],
            key=lambda op: op.get('cond', 0)
        )

        if not sorted_ops:
            lines.append("Keine Konditionszahlen verfügbar.")
            report = "\n".join(lines)
            self._logger.info(report)
            return report

        # Jede Operation auflisten
        for op in sorted_ops:
            name = op.get('name', 'Unbekannt')
            cond = op.get('cond', float('nan'))
            operation = op.get('operation', '')

            # Bewertung nach Konditionszahl
            if cond < 1e3:
                status = "gut konditioniert"
            elif cond < 1e8:
                status = "mäßig konditioniert"
            elif cond < 1e12:
                status = "schlecht konditioniert"
            else:
                status = "KRITISCH (fast singulär)"

            entry = f"  {name} [{operation}]: κ = {cond:.2e} → {status}"
            lines.append(entry)

            # Auch im Logger ausgeben
            if cond > 1e10:
                self._logger.warning(entry.strip())
            else:
                self._logger.info(entry.strip())

        lines.append("=========================")
        report = "\n".join(lines)
        return report

    def disable(self) -> None:
        """
        @brief Deaktiviert den Logger komplett (kein Output mehr).
        @description
            Setzt den Level auf CRITICAL+1, sodass keine Meldung mehr ausgegeben wird.
            Nützlich für Produktionsbetrieb oder Performance-Messungen.

        @author Michael Fuhrmann
        @lastModified 2026-03-10
        """
        self._logger.setLevel(logging.CRITICAL + 1)


# ===========================================================================
# GLOBALER STANDARD-LOGGER UND HILFSFUNKTIONEN
# ===========================================================================

# Globaler Singleton-Logger für projektweite Nutzung ohne Instanziierung
_default_logger = MathLogger(name="default", level="INFO")


def get_logger(name: str = None) -> MathLogger:
    """
    @brief Gibt den globalen Logger zurück oder erstellt einen neuen benannten Logger.
    @description
        Ohne Argument: Gibt den globalen Standard-Logger zurück.
        Mit Name: Erstellt einen neuen Logger mit diesem Namen.

    @param name: Optionaler Name für einen eigenen Logger
    @return MathLogger-Instanz
    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    global _default_logger
    if name is None:
        return _default_logger
    # Neuen Logger mit gegebenem Namen erstellen
    return MathLogger(name=name)


def enable_debug_logging() -> None:
    """
    @brief Aktiviert DEBUG-Level für den globalen Logger.
    @description
        Schaltet alle Berechnungsschritte auf Konsole ein.
        Nützlich für die Entwicklung und Fehlersuche.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    global _default_logger
    _default_logger.set_level("DEBUG")


def disable_logging() -> None:
    """
    @brief Deaktiviert Logging des globalen Loggers komplett.
    @description
        Unterdrückt alle Ausgaben für Performance-Messungen oder saubere Ausgabe.

    @author Michael Fuhrmann
    @lastModified 2026-03-10
    """
    global _default_logger
    _default_logger.disable()
