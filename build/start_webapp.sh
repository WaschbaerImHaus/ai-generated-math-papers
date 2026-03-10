#!/bin/bash
#
# @file start_webapp.sh
# @brief Startet die Flask-Webanwendung des Mathematik-Spezialisten.
# @description
#     Wechselt in das webapp-Verzeichnis und startet die Flask-App
#     auf Port 8080. Die PID wird in einer Datei gespeichert,
#     damit stop_webapp.sh den Prozess beenden kann.
# @author Kurt Ingwer
# @date 2026-03-10
#

# Verzeichnis dieser Datei ermitteln (unabhängig vom Aufrufpfad)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WEBAPP_DIR="$SCRIPT_DIR/../src/webapp"
PID_FILE="$SCRIPT_DIR/specialist-maths.pid"

# Prüfen ob die Webapp schon läuft
if [ -f "$PID_FILE" ]; then
    OLD_PID=$(cat "$PID_FILE")
    if kill -0 "$OLD_PID" 2>/dev/null; then
        echo "[specialist-maths] Webapp läuft bereits (PID $OLD_PID)"
        echo "  → Zum Neustarten: $SCRIPT_DIR/stop_webapp.sh && $SCRIPT_DIR/start_webapp.sh"
        exit 1
    else
        # PID-Datei veraltet (Prozess nicht mehr aktiv)
        rm -f "$PID_FILE"
    fi
fi

# Webapp starten
echo "====================================="
echo "  Mathematik-Spezialist Web-Interface"
echo "  Build 12 | Port 8080"
echo "  URL: http://localhost:8080"
echo "====================================="

cd "$WEBAPP_DIR" || { echo "Fehler: Verzeichnis $WEBAPP_DIR nicht gefunden."; exit 1; }

# Im Hintergrund starten, PID speichern
python3 app.py &
APP_PID=$!
echo $APP_PID > "$PID_FILE"

echo "[specialist-maths] Gestartet (PID $APP_PID)"
echo "[specialist-maths] PID gespeichert: $PID_FILE"
echo "[specialist-maths] Zum Stoppen: $SCRIPT_DIR/stop_webapp.sh"
