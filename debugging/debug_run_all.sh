#!/bin/bash
# @file debug_run_all.sh
# @brief Führt alle Debug-Skripte nacheinander aus.
# @author Kurt Ingwer
# @date 2026-03-07

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "========================================"
echo "specialist-maths - Debug-Suite"
echo "========================================"

for script in "$SCRIPT_DIR"/debug_*.py; do
    echo ""
    echo ">>> $(basename "$script")"
    python3 "$script"
done

echo ""
echo "========================================"
echo "Alle Debug-Skripte abgeschlossen."
