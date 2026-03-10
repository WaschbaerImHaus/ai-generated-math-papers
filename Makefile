# Makefile — specialist-maths
# Zentrale Testkonfiguration: pytest + Hypothesis + Doctests in einem Lauf
# Autor: Kurt Ingwer
# Zuletzt geändert: 2026-03-10

.PHONY: test test-fast test-doctests test-property test-coverage clean lint

# --- Vollständiger Testlauf (Standard) ---
# Führt alle Tests aus: Unit, Hypothesis-Property-Tests, Doctests
test:
	@echo "=== Vollständiger Testlauf ==="
	cd $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) && \
	python3 -m pytest tests/ \
	    --doctest-modules \
	    --ignore=src/webapp \
	    --ignore=src/__init__.py \
	    -p hypothesis \
	    --hypothesis-seed=0 \
	    -v \
	    --tb=short \
	    -q \
	    2>&1

# --- Schneller Testlauf (kein Doctest, kein Hypothesis) ---
test-fast:
	@echo "=== Schneller Testlauf (ohne Doctests & Property-Tests) ==="
	cd $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) && \
	python3 -m pytest tests/ \
	    -p no:hypothesis \
	    --tb=short \
	    -q \
	    2>&1

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
	    -v \
	    --tb=short \
	    2>&1

# --- Testabdeckung (Coverage) ---
test-coverage:
	@echo "=== Test-Coverage ==="
	cd $(dir $(abspath $(lastword $(MAKEFILE_LIST)))) && \
	python3 -m pytest tests/ \
	    --cov=src \
	    --cov-report=term-missing \
	    --cov-report=html:build/coverage \
	    --ignore=src/webapp \
	    --tb=short \
	    -q \
	    2>&1
	@echo "Coverage-Bericht: build/coverage/index.html"

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
