#!/bin/bash
# =============================================================================
# @file run_doctests.sh
# @brief Führt alle Doctests in src/ aus und zeigt eine Zusammenfassung.
# @description
#     Durchläuft alle Python-Dateien im src/-Verzeichnis und führt für jede
#     Datei "python3 -m doctest" aus. Am Ende wird eine Gesamtübersicht
#     mit Erfolg/Fehlschlag-Statistik ausgegeben.
#
#     Verwendung:
#         bash debugging/run_doctests.sh
#     oder direkt:
#         ./debugging/run_doctests.sh
#
#     Exitcode:
#     - 0: Alle Doctests bestanden
#     - 1: Mindestens ein Doctest fehlgeschlagen
#
# @author Kurt Ingwer
# @lastModified 2026-03-10
# =============================================================================

# Skript bei Fehlern NICHT automatisch abbrechen – wir wollen alle Dateien prüfen
set +e

# Absoluter Pfad zum src/-Verzeichnis (relativ zur Position dieses Skripts)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_DIR="$SCRIPT_DIR/../src"

# Farben für die Ausgabe (nur wenn Terminal Farben unterstützt)
if [ -t 1 ]; then
    GREEN="\033[0;32m"
    RED="\033[0;31m"
    YELLOW="\033[0;33m"
    RESET="\033[0m"
else
    GREEN=""
    RED=""
    YELLOW=""
    RESET=""
fi

echo "=== Doctests für specialist-maths ==="
echo "Verzeichnis: $SRC_DIR"
echo "Datum: $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================="

# Zähler für Statistik
TOTAL=0
PASSED=0
FAILED=0

# Alle Python-Dateien im src/-Verzeichnis durchlaufen
for f in "$SRC_DIR"/*.py; do
    # Nur reguläre Dateien verarbeiten
    [ -f "$f" ] || continue

    # Dateiname ohne Pfad für Ausgabe
    basename_f="$(basename "$f")"

    echo ""
    echo "--- $basename_f ---"

    # Doctest ausführen und letzten 5 Zeilen der Ausgabe anzeigen
    # -v für verbose zeigt jedes einzelne Beispiel
    OUTPUT=$(cd "$SRC_DIR" && python3 -m doctest "$basename_f" -v 2>&1)
    EXIT_CODE=$?

    # Nur die letzten 5 Zeilen der Ausgabe anzeigen (Zusammenfassung)
    echo "$OUTPUT" | tail -5

    TOTAL=$((TOTAL + 1))

    if [ $EXIT_CODE -eq 0 ]; then
        echo -e "${GREEN}[OK] $basename_f${RESET}"
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}[FEHLER] $basename_f – Exit-Code: $EXIT_CODE${RESET}"
        FAILED=$((FAILED + 1))
    fi
done

# Abschlussstatistik ausgeben
echo ""
echo "======================================="
echo "=== Zusammenfassung ==="
echo "Gesamt:      $TOTAL Dateien"
echo -e "Bestanden:   ${GREEN}$PASSED${RESET}"
echo -e "Fehlgeschlagen: ${RED}$FAILED${RESET}"
echo "======================================="

# Exitcode: 0 wenn alle erfolgreich, sonst 1
if [ $FAILED -gt 0 ]; then
    echo -e "${RED}ERGEBNIS: Doctests fehlgeschlagen!${RESET}"
    exit 1
else
    echo -e "${GREEN}ERGEBNIS: Alle Doctests bestanden.${RESET}"
    exit 0
fi
