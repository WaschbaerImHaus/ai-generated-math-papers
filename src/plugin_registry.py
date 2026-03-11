"""
@file plugin_registry.py
@brief Plugin-Registrierungs-System für mathematische Module.
@description
    Ermöglicht die dynamische Registrierung, das Laden und die Verwaltung
    von Plugin-Modulen. Plugins werden als Name → Modulpfad-Paare registriert
    und bei Bedarf via importlib.import_module geladen.

    Verwendung:
        from plugin_registry import DEFAULT_REGISTRY

        # Plugin registrieren
        DEFAULT_REGISTRY.register("mein_modul", "mein_paket.mein_modul")

        # Plugin laden
        modul = DEFAULT_REGISTRY.load("mein_modul")

        # Alle Plugins auflisten
        print(DEFAULT_REGISTRY.list_plugins())

    Das System enthält eine vorbefüllte DEFAULT_REGISTRY mit allen
    Standardmodulen des specialist-maths-Projekts.

@author Michael Fuhrmann
@date 2026-03-11
@lastModified 2026-03-11
"""

import importlib
import importlib.util
import sys
import os
from typing import Any, Optional


# ---------------------------------------------------------------------------
# KLASSE: PluginRegistry
# ---------------------------------------------------------------------------

class PluginRegistry:
    """
    @brief Registrierung und dynamisches Laden von Plugin-Modulen.
    @description
        Verwaltet eine Tabelle aus Plugin-Namen und Modulpfaden.
        Plugins können zur Laufzeit registriert, geladen und entfernt werden.

        Intern wird ein Dictionary {name: modulpfad} für die Registrierung
        und ein weiteres {name: modul} für den Loaded-Cache verwendet.

    @author Michael Fuhrmann
    @lastModified 2026-03-11
    """

    def __init__(self) -> None:
        """
        @brief Initialisiert eine leere Plugin-Registrierung.
        @description
            Erstellt zwei interne Dicts:
            - _registry: Speichert Name → Modulpfad-Zuordnungen
            - _loaded:   Speichert bereits geladene Modul-Objekte (Cache)

        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Registrierungstabelle: Plugin-Name → Modulpfad
        self._registry: dict[str, str] = {}

        # Modul-Cache: Plugin-Name → bereits importiertes Modul-Objekt
        self._loaded: dict[str, Any] = {}

    def register(self, name: str, module_path: str) -> None:
        """
        @brief Registriert ein Plugin unter dem angegebenen Namen.
        @description
            Speichert die Zuordnung Name → Modulpfad in der internen
            Registrierungstabelle. Falls ein Plugin mit diesem Namen
            bereits registriert ist, wird der Modulpfad überschrieben
            und der eventuell gecachte Import wird verworfen.

        @param name: Eindeutiger Name des Plugins (z.B. "algebra_core").
        @param module_path: Importpfad des Moduls (z.B. "algebra_numbertheory").
        @raises ValueError: Wenn Name oder Modulpfad leer ist.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Eingabevalidierung: Name und Pfad dürfen nicht leer sein
        if not name or not name.strip():
            raise ValueError("Plugin-Name darf nicht leer sein")
        if not module_path or not module_path.strip():
            raise ValueError(f"Modulpfad für Plugin '{name}' darf nicht leer sein")

        # Plugin registrieren (überschreibt ggf. bestehenden Eintrag)
        self._registry[name] = module_path

        # Gecachtes Modul verwerfen (Neuladen erzwingen)
        if name in self._loaded:
            del self._loaded[name]

    def load(self, name: str) -> Any:
        """
        @brief Lädt ein registriertes Plugin-Modul und gibt es zurück.
        @description
            Falls das Plugin bereits geladen wurde, wird das gecachte
            Modul-Objekt zurückgegeben (kein erneuter Import).
            Andernfalls wird das Modul via importlib.import_module geladen
            und im Cache gespeichert.

            Suchpfad: Das src/-Verzeichnis wird automatisch zu sys.path
            hinzugefügt, damit alle Projektmodule gefunden werden.

        @param name: Name des zu ladenden Plugins.
        @return: Das geladene Modul-Objekt.
        @raises KeyError: Wenn das Plugin nicht registriert ist.
        @raises ImportError: Wenn das Modul nicht importiert werden kann.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Prüfen ob Plugin registriert ist
        if name not in self._registry:
            raise KeyError(
                f"Plugin '{name}' ist nicht registriert. "
                f"Verfügbare Plugins: {list(self._registry.keys())}"
            )

        # Gecachtes Modul zurückgeben falls schon geladen
        if name in self._loaded:
            return self._loaded[name]

        # Modulpfad aus der Registrierung holen
        module_path = self._registry[name]

        # src/-Verzeichnis zum Suchpfad hinzufügen (für Projektmodule)
        src_dir = os.path.dirname(os.path.abspath(__file__))
        if src_dir not in sys.path:
            sys.path.insert(0, src_dir)

        # Modul via importlib laden
        module = importlib.import_module(module_path)

        # Im Cache speichern für schnellere Folgeaufrufe
        self._loaded[name] = module
        return module

    def list_plugins(self) -> list[str]:
        """
        @brief Gibt eine Liste aller registrierten Plugin-Namen zurück.
        @description
            Die Liste ist alphabetisch sortiert für konsistente Ausgabe.
            Gibt alle registrierten Plugins zurück – unabhängig davon,
            ob sie bereits geladen wurden.

        @return: Alphabetisch sortierte Liste der Plugin-Namen.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Alphabetisch sortieren für konsistente Reihenfolge
        return sorted(self._registry.keys())

    def is_loaded(self, name: str) -> bool:
        """
        @brief Prüft ob ein Plugin bereits geladen (importiert) wurde.
        @description
            Gibt True zurück, wenn das Plugin mindestens einmal via load()
            erfolgreich geladen wurde und sich im Modul-Cache befindet.
            Gibt False zurück, wenn das Plugin (noch) nicht geladen wurde
            oder nicht registriert ist.

        @param name: Name des Plugins.
        @return: True wenn geladen, sonst False.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Im Loaded-Cache nachschauen
        return name in self._loaded

    def unregister(self, name: str) -> None:
        """
        @brief Entfernt ein Plugin aus der Registrierung.
        @description
            Löscht sowohl den Registrierungs-Eintrag als auch das
            eventuell gecachte Modul-Objekt. Nach dem Aufruf ist das
            Plugin nicht mehr über dieses Registry zugänglich.

            Falls das Plugin nicht registriert ist, wird kein Fehler
            geworfen (idempotente Operation).

        @param name: Name des zu entfernenden Plugins.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Registrierungs-Eintrag entfernen (kein Fehler wenn nicht vorhanden)
        self._registry.pop(name, None)

        # Gecachtes Modul-Objekt entfernen
        self._loaded.pop(name, None)

    def reload(self, name: str) -> Any:
        """
        @brief Lädt ein Plugin-Modul neu (erzwingt erneuten Import).
        @description
            Löscht den Modul-Cache für das angegebene Plugin und
            importiert es erneut. Nützlich wenn sich das Modul geändert hat.

        @param name: Name des neu zu ladenden Plugins.
        @return: Das neu geladene Modul-Objekt.
        @raises KeyError: Wenn das Plugin nicht registriert ist.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        # Cache-Eintrag löschen (erzwingt Neuladen)
        self._loaded.pop(name, None)

        # Neu laden via load()
        return self.load(name)

    def __repr__(self) -> str:
        """
        @brief Lesbare String-Darstellung der Registry.
        @return: String mit Anzahl registrierter und geladener Plugins.
        @author Michael Fuhrmann
        @lastModified 2026-03-11
        """
        return (
            f"PluginRegistry("
            f"registriert={len(self._registry)}, "
            f"geladen={len(self._loaded)})"
        )


# ---------------------------------------------------------------------------
# GLOBALE STANDARD-REGISTRY mit vorbefüllten Projektmodulen
# ---------------------------------------------------------------------------

#: Globale Standard-Registry für das specialist-maths-Projekt.
#: Enthält alle Kern-Module des Projekts als vorregistrierte Plugins.
DEFAULT_REGISTRY = PluginRegistry()

# Alle Kern-Module des Projekts vorregistrieren
# Schlüssel: Plugin-Name (kurz, einprägsam)
# Wert: Python-Importpfad (wie in import-Anweisung)

# --- Algebra und Zahlentheorie ---
DEFAULT_REGISTRY.register("algebra_core",           "algebra")
DEFAULT_REGISTRY.register("algebra_numbertheory",   "algebra_numbertheory")

# --- Analysis ---
DEFAULT_REGISTRY.register("analysis",               "analysis")

# --- Lineare Algebra ---
DEFAULT_REGISTRY.register("linear_algebra",         "linear_algebra")
DEFAULT_REGISTRY.register("matrix_ops",             "matrix_ops")

# --- Fourier-Analyse ---
DEFAULT_REGISTRY.register("fourier",                "fourier")

# --- Statistik und Stochastik ---
DEFAULT_REGISTRY.register("statistics",             "statistics_math")

# --- Differentialgleichungen ---
DEFAULT_REGISTRY.register("ode",                    "ode")

# --- Zahlentheorie (erweitert) ---
DEFAULT_REGISTRY.register("proof_theory",           "proof_theory")
DEFAULT_REGISTRY.register("analytic_number_theory", "analytic_number_theory")

# --- Tensor- und Differentialgeometrie ---
DEFAULT_REGISTRY.register("tensor_geometry",        "tensor_geometry")

# --- Visualisierung ---
DEFAULT_REGISTRY.register("visualization",          "visualization")

# --- SageMath-Interoperabilität ---
DEFAULT_REGISTRY.register("sagemath_bridge",        "sagemath_bridge")

# --- REPL-Interface ---
DEFAULT_REGISTRY.register("repl",                   "repl")

# --- Kategorientheorie ---
DEFAULT_REGISTRY.register("category_theory",        "category_theory")

# --- Algebraische Strukturen ---
DEFAULT_REGISTRY.register("algebraic_structures",   "algebraic_structures")
DEFAULT_REGISTRY.register("algebraic_topology",     "algebraic_topology")
DEFAULT_REGISTRY.register("algebraic_geometry",     "algebraic_geometry")

# --- Galois & Kommutative Algebra ---
DEFAULT_REGISTRY.register("galois_theory",          "galois_theory")
DEFAULT_REGISTRY.register("commutative_algebra",    "commutative_algebra")
DEFAULT_REGISTRY.register("ring_theory",            "ring_theory")
DEFAULT_REGISTRY.register("group_theory",           "group_theory")

# --- Zahlentheorie (erweitert) ---
DEFAULT_REGISTRY.register("p_adic",                 "p_adic")
DEFAULT_REGISTRY.register("modular_forms",          "modular_forms")
DEFAULT_REGISTRY.register("l_functions",            "l_functions")
DEFAULT_REGISTRY.register("iwasawa_theory",         "iwasawa_theory")

# --- Geometrie ---
DEFAULT_REGISTRY.register("differential_geometry",  "differential_geometry")
DEFAULT_REGISTRY.register("symplectic_geometry",    "symplectic_geometry")
DEFAULT_REGISTRY.register("classical_geometry",     "classical_geometry")

# --- Logik & Grundlagen ---
DEFAULT_REGISTRY.register("mathematical_logic",     "mathematical_logic")
DEFAULT_REGISTRY.register("set_theory",             "set_theory")
DEFAULT_REGISTRY.register("model_theory",           "model_theory")

# --- Funktionalanalysis & Maßtheorie ---
DEFAULT_REGISTRY.register("functional_analysis",    "functional_analysis")
DEFAULT_REGISTRY.register("measure_theory",         "measure_theory")
DEFAULT_REGISTRY.register("operator_algebras",      "operator_algebras")

# --- Kombinatorik & Diskrete Mathematik ---
DEFAULT_REGISTRY.register("combinatorics",          "combinatorics")
DEFAULT_REGISTRY.register("graph_theory",           "graph_theory")
DEFAULT_REGISTRY.register("coding_theory",          "coding_theory")

# --- Höhere Präzision ---
DEFAULT_REGISTRY.register("arbitrary_precision",    "arbitrary_precision")
