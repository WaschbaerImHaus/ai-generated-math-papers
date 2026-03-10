@echo off
REM start_server.bat — Startet den Mathematik-Spezialisten Web-Server
REM Autor: Kurt Ingwer
REM Zuletzt geändert: 2026-03-10

title Mathematik-Spezialist Web-Server (Port 8080)

REM Ins Projektverzeichnis wechseln (Skript liegt im Projektstamm)
cd /d "%~dp0"

REM Python prüfen
python --version >nul 2>&1
if errorlevel 1 (
    echo FEHLER: Python nicht gefunden. Bitte Python 3.10+ installieren.
    pause
    exit /b 1
)

REM Flask prüfen
python -c "import flask" >nul 2>&1
if errorlevel 1 (
    echo Flask nicht gefunden. Installiere Abhängigkeiten...
    pip install flask
)

echo.
echo  =============================================
echo   Mathematik-Spezialist Web-Server
echo   http://localhost:8080
echo   Strg+C zum Beenden
echo  =============================================
echo.

REM Server starten
set PYTHONPATH=%~dp0src
python src\webapp\app.py

pause
