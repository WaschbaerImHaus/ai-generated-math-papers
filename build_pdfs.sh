#!/bin/bash
# PDF-Kompilierung aller LaTeX-Papers
# Autor: Michael Fuhrmann
# Erzeugt PDFs unter papers-pdf/ (nicht im Git-Repo)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPERS_DIR="$SCRIPT_DIR/papers"
OUT_DIR="$SCRIPT_DIR/papers-pdf"
LOG_FILE="$OUT_DIR/build.log"

# Timeout pro pdflatex-Aufruf in Sekunden (60s verhindert hängende TikZ/etc.)
PDFLATEX_TIMEOUT=60

echo "=== specialist-maths PDF Build ==="
echo "Quelle: $PAPERS_DIR"
echo "Ausgabe: $OUT_DIR"
echo ""

mkdir -p "$OUT_DIR"
echo "" > "$LOG_FILE"

SUCCESS=0
FAILED=0
FAILED_FILES=()

# Alle .tex-Dateien finden
while IFS= read -r -d '' TEX_FILE; do
    BASENAME=$(basename "$TEX_FILE" .tex)
    # Relativer Pfad ab papers/
    REL=$(realpath --relative-to="$PAPERS_DIR" "$(dirname "$TEX_FILE")")
    TARGET_DIR="$OUT_DIR/$REL"
    mkdir -p "$TARGET_DIR"

    echo -n "  Kompiliere: $REL/$BASENAME.tex ... "

    # Zweimal kompilieren für korrekte Referenzen (mit Timeout gegen hängende Prozesse)
    if timeout "$PDFLATEX_TIMEOUT" pdflatex -interaction=nonstopmode -halt-on-error \
        -output-directory="$TARGET_DIR" \
        "$TEX_FILE" >> "$LOG_FILE" 2>&1 && \
       timeout "$PDFLATEX_TIMEOUT" pdflatex -interaction=nonstopmode -halt-on-error \
        -output-directory="$TARGET_DIR" \
        "$TEX_FILE" >> "$LOG_FILE" 2>&1; then
        echo "OK"
        SUCCESS=$((SUCCESS + 1))
        # Hilfsdateien aufräumen, nur PDF behalten
        rm -f "$TARGET_DIR/$BASENAME.aux" \
              "$TARGET_DIR/$BASENAME.log" \
              "$TARGET_DIR/$BASENAME.out" \
              "$TARGET_DIR/$BASENAME.toc" \
              "$TARGET_DIR/$BASENAME.fls" \
              "$TARGET_DIR/$BASENAME.fdb_latexmk"
    else
        echo "FEHLER"
        FAILED=$((FAILED + 1))
        FAILED_FILES+=("$REL/$BASENAME.tex")
    fi
done < <(find "$PAPERS_DIR" -name "*.tex" -print0 | sort -z)

echo ""
echo "=== Ergebnis ==="
echo "  Erfolgreich: $SUCCESS"
echo "  Fehler:      $FAILED"

if [ ${#FAILED_FILES[@]} -gt 0 ]; then
    echo ""
    echo "  Fehlgeschlagene Dateien:"
    for f in "${FAILED_FILES[@]}"; do
        echo "    - $f"
    done
fi

echo ""
echo "PDFs liegen in: $OUT_DIR"
echo "Build-Log:      $LOG_FILE"
