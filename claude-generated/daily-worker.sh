#!/bin/bash
##
# @file daily-worker.sh
# @brief Täglicher automatischer Arbeits-Lauf für das Projekt specialist-maths.
#
# @description
#   Wird täglich um 07:00 Uhr per Cronjob ausgeführt.
#   Startet Claude Code im Projektverzeichnis mit dem Auftrag,
#   offene Aufgaben gemäß CLAUDE.md fortzusetzen:
#     - MEMORY.md, FEATURES.md, BUGS.md lesen
#     - nächste geplante Features implementieren (TDD)
#     - Tests ausführen
#     - Build hochzählen, committen und pushen
#
#   Verwendet flock, damit nie zwei Instanzen gleichzeitig laufen.
#   Log-Ausgabe in /home/claude-code/.claude/specialist-maths-daily.log
#
# @author Kurt Ingwer
# @date   2026-03-08
# @lastModified: 2026-03-08
##

# --- Root-Schutz --------------------------------------------------------------
if [ "$(id -u)" -eq 0 ]; then
    echo "FEHLER: Dieses Skript darf nicht als root ausgeführt werden." >&2
    exit 1
fi

# --- Konfiguration ------------------------------------------------------------
PROJECT_DIR="/home/claude-code/project/specialist-maths"
CLAUDE_BIN="/home/claude-code/.local/bin/claude"
LOGFILE="/home/claude-code/.claude/specialist-maths-daily.log"
LOCKFILE="/tmp/specialist-maths-daily.lock"

# Zeitstempel für Log
echo "========================================" >> "$LOGFILE"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Start daily-worker specialist-maths" >> "$LOGFILE"

# --- Lockfile (verhindert parallele Ausführung) -------------------------------
exec 9>"$LOCKFILE"
if ! flock -n 9; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Bereits aktiv (Lock vorhanden). Abbruch." >> "$LOGFILE"
    exit 0
fi

# --- Voraussetzungen prüfen ---------------------------------------------------
if [ ! -f "$CLAUDE_BIN" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - FEHLER: Claude Binary nicht gefunden: $CLAUDE_BIN" >> "$LOGFILE"
    exit 1
fi

if [ ! -d "$PROJECT_DIR" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - FEHLER: Projektverzeichnis nicht gefunden: $PROJECT_DIR" >> "$LOGFILE"
    exit 1
fi

# --- Claude ausführen ---------------------------------------------------------
cd "$PROJECT_DIR" || exit 1

"$CLAUDE_BIN" \
    --dangerously-skip-permissions \
    --continue \
    -p "Handle nach der [home]/.claude/CLAUDE.md. erledige unerledigtes für das Projekt specialist-maths." \
    >> "$LOGFILE" 2>&1

EXIT_CODE=$?

echo "$(date '+%Y-%m-%d %H:%M:%S') - Beendet mit Exit-Code $EXIT_CODE" >> "$LOGFILE"
echo "========================================" >> "$LOGFILE"

exit $EXIT_CODE
