#!/bin/bash
#
# @file stop_webapp.sh
# @brief Beendet die Flask-Webanwendung des Mathematik-Spezialisten.
# @description
#     Liest die gespeicherte PID aus der .pid-Datei und beendet
#     den Flask-Prozess. Falls keine PID-Datei vorhanden ist,
#     wird nach dem Prozess per Name gesucht (Fallback).
# @author Kurt Ingwer
# @date 2026-03-10
#

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PID_FILE="$SCRIPT_DIR/specialist-maths.pid"

if [ -f "$PID_FILE" ]; then
    # PID-Datei vorhanden → Prozess direkt beenden
    PID=$(cat "$PID_FILE")
    if kill -0 "$PID" 2>/dev/null; then
        echo "[specialist-maths] Beende Prozess (PID $PID)..."
        kill "$PID"
        sleep 1
        # Falls noch nicht beendet, SIGKILL schicken
        if kill -0 "$PID" 2>/dev/null; then
            echo "[specialist-maths] Erzwinge Beenden (SIGKILL)..."
            kill -9 "$PID"
        fi
        rm -f "$PID_FILE"
        echo "[specialist-maths] Gestoppt."
    else
        echo "[specialist-maths] Prozess $PID ist nicht mehr aktiv."
        rm -f "$PID_FILE"
    fi
else
    # Fallback: Prozess per Name suchen
    echo "[specialist-maths] Keine PID-Datei gefunden. Suche per Prozessname..."
    PIDS=$(pgrep -f "python3.*app.py" 2>/dev/null)
    if [ -n "$PIDS" ]; then
        echo "[specialist-maths] Gefundene Prozesse: $PIDS"
        echo "$PIDS" | xargs kill 2>/dev/null
        echo "[specialist-maths] Beendet."
    else
        echo "[specialist-maths] Kein laufender Prozess gefunden."
    fi
fi
