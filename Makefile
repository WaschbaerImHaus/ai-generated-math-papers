# Makefile — specialist-maths
# Zentrale Testkonfiguration: pytest + xdist (parallel) + Coverage
# Autor: Kurt Ingwer
# Zuletzt geändert: 2026-03-11

.PHONY: test test-fast test-doctests test-property test-coverage clean lint

# --- Vollständiger Testlauf mit paralleler Ausführung (Standard) ---
# -n auto: Alle verfügbaren CPU-Kerne nutzen (pytest-xdist)
# --timeout=60: Einzelner Test darf maximal 60 Sekunden dauern
test:
	cd /home/claude-code/project/specialist-maths && python3 -m pytest tests/ -n auto --timeout=60

# --- Schneller Testlauf: erster Fehler stoppt (-x), minimale Ausgabe (-q) ---
test-fast:
	cd /home/claude-code/project/specialist-maths && python3 -m pytest tests/ -n auto --timeout=60 -x -q

# --- Testabdeckung (Coverage) ---
# Erstellt HTML-Bericht unter build/coverage/index.html
test-coverage:
	cd /home/claude-code/project/specialist-maths && python3 -m pytest tests/ --cov=src --cov-report=html

# --- Nur Doctests ---
test-doctests:
	@echo "=== Doctest-Lauf ==="
	cd $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) && \
	python3 -m pytest src/ \
	    --doctest-modules \
	    --ignore=src/webapp \
	    --ignore=src/__init__.py \
	    --tb=short \
	    -q \
	    2>&1

# --- Nur Hypothesis Property-Based Tests ---
test-property:
	@echo "=== Property-Based Tests (Hypothesis) ==="
	cd $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) && \
	python3 -m pytest tests/ \
	    -k "property or hypothesis or fuzz" \
	    --hypothesis-seed=0 \
	    -n auto \
	    -v \
	    --tb=short \
	    2>&1

# --- Lint (SyntaxWarnings prüfen) ---
lint:
	@echo "=== SyntaxWarning-Prüfung ==="
	cd $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) && \
	python3 -W error -m py_compile src/*.py && echo "Keine SyntaxWarnings gefunden."

# --- Build-Infos ---
status:
	@echo "Build: $$(cat src/build.txt)"
	@echo "Tests: $$(python3 -m pytest tests/ --collect-only -q 2>&1 | tail -1)"
	@echo "Python: $$(python3 --version)"

# --- Aufräumen ---
clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
	rm -rf .pytest_cache .hypothesis build/coverage 2>/dev/null || true
	@echo "Aufgeräumt."
