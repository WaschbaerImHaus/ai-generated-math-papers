"""
@file __init__.py
@brief Lazy-Loading Package für specialist-maths.
@description
    Module werden erst beim ersten Zugriff geladen (Lazy-Loading via importlib).
    Das reduziert die Startzeit erheblich, wenn nur einzelne Module benötigt werden.

    Verwendung (gleich wie vorher):
        from src import algebra
        from src.algebra import gcd       # lädt nur algebra
        from src.linear_algebra import Vector

    Nicht geladene Module belegen keinen Speicher und verursachen
    keine Import-Nebeneffekte, bis sie tatsächlich benötigt werden.

    Verfügbare Module können mit list_modules() abgerufen werden.

@author Kurt Ingwer
@date 2026-03-10
@lastModified 2026-03-10
"""

import importlib
import sys
from typing import Any

# ---------------------------------------------------------------------------
# Verzeichnis der verfügbaren Module
# Schlüssel: Attributname, Wert: Modulname (relativ zum Package)
# ---------------------------------------------------------------------------
_MODULE_MAP: dict[str, str] = {
    # Algebra und Zahlentheorie
    'algebra':               'algebra',
    'algebra_core':          'algebra_core',
    'algebra_numbertheory':  'algebra_numbertheory',
    'algebra_diophantine':   'algebra_diophantine',

    # Analysis und Differentialrechnung
    'analysis':              'analysis',

    # Lineare Algebra – aufgeteilt in Teilmodule
    'linear_algebra':        'linear_algebra',
    'vectors':               'vectors',
    'matrix_ops':            'matrix_ops',
    'matrix_decomp':         'matrix_decomp',

    # Statistik und Wahrscheinlichkeitsrechnung
    'statistics_math':       'statistics_math',

    # Differentialgleichungen
    'ode':                   'ode',

    # Fourier-Analysis und Spektrum
    'fourier':               'fourier',

    # Komplexe Analysis und Riemannsche Flächen
    'complex_analysis':      'complex_analysis',

    # Numerische Methoden – Interpolation, Optimierung
    'numerical_methods':     'numerical_methods',

    # Modulformen und elliptische Kurven
    'modular_forms':         'modular_forms',

    # P-adische Zahlen und nicht-archimedische Normen
    'p_adic':                'p_adic',

    # Visualisierung (Matplotlib-basiert)
    'visualization':         'visualization',

    # Beweistheorie und Zahlentheoretische Vermutungen
    'proof_theory':          'proof_theory',

    # Analytische Zahlentheorie und Primzahlverteilung
    'analytic_number_theory': 'analytic_number_theory',

    # Topologie, metrische Räume, Mannigfaltigkeiten
    'topology':              'topology',

    # Graphentheorie und Kombinatorik
    'graph_theory':          'graph_theory',

    # Millennium-Probleme (experimentell)
    'millennium_problems':   'millennium_problems',

    # LaTeX-Export mathematischer Ausdrücke
    'latex_export':          'latex_export',

    # REPL-Modus (interaktive Eingabe)
    'repl':                  'repl',

    # Konfiguration und Einstellungen
    'config':                'config',

    # Ausnahmen und Fehlerklassen
    'exceptions':            'exceptions',

    # Logging-Infrastruktur
    'math_logger':           'math_logger',

    # Schrittweise Berechnungsanzeige
    'step_by_step':          'step_by_step',
}


def __getattr__(name: str) -> Any:
    """
    @brief Lazy-Import: Modul wird erst beim ersten Zugriff geladen.
    @description
        Python ruft diese Funktion auf, wenn ein Attribut nicht direkt
        im Package-Namespace gefunden wird. Wir laden dann das Modul
        via importlib und cachen es im globalen Namespace.

        Mechanismus:
        1. Zugriff auf 'specialist_maths.algebra'
        2. Python findet 'algebra' nicht in globals()
        3. __getattr__('algebra') wird aufgerufen
        4. importlib.import_module('.algebra', package) lädt das Modul
        5. globals()['algebra'] = module (für künftige Zugriffe gecacht)
        6. Das Modul wird zurückgegeben

    @param name Attributname, der aufgerufen wurde.
    @return Geladenes Modul.
    @raises AttributeError Wenn das Modul nicht im _MODULE_MAP-Verzeichnis steht.
    @lastModified 2026-03-10
    """
    if name in _MODULE_MAP:
        # Modulname aus der Karte lesen
        module_name = _MODULE_MAP[name]
        # Relativer Import innerhalb des src-Packages
        module = importlib.import_module(f'.{module_name}', package=__name__)
        # Im globalen Namespace cachen, damit spätere Zugriffe direkt gehen
        globals()[name] = module
        return module

    # Unbekanntes Attribut: aussagekräftige Fehlermeldung mit allen Optionen
    raise AttributeError(
        f"Modul '{name}' nicht in specialist-maths gefunden. "
        f"Verfügbare Module: {sorted(_MODULE_MAP.keys())}"
    )


def list_modules() -> list[str]:
    """
    @brief Gibt eine sortierte Liste aller verfügbaren Module zurück.
    @description
        Nützlich zur Übersicht und für dynamische Importe.
        Gibt alle registrierten Module zurück, unabhängig davon ob sie
        bereits geladen wurden.

    @return Sortierte Liste der Modulnamen.
    @lastModified 2026-03-10
    """
    return sorted(_MODULE_MAP.keys())


def loaded_modules() -> list[str]:
    """
    @brief Gibt eine Liste der bereits geladenen Module zurück.
    @description
        Zeigt welche Module tatsächlich importiert wurden.
        Nützlich zur Diagnose und Performance-Analyse.

    @return Liste der Namen bereits geladener Module.
    @lastModified 2026-03-10
    """
    # Prüft welche Modulnamen bereits im globals()-Dict stehen
    return [name for name in _MODULE_MAP if name in globals()]


# ---------------------------------------------------------------------------
# Versionsinformation
# ---------------------------------------------------------------------------
__version__: str = "13.0"
__author__: str = "Kurt Ingwer"
__all__: list[str] = list(_MODULE_MAP.keys()) + ['list_modules', 'loaded_modules']
