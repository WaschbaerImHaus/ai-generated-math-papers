# SECURITY_RISKS.md - specialist-maths

## Offene Sicherheitsrisiken

### MEDIUM: sp.sympify() mit Nutzereingaben (analysis.py)
**Entdeckt:** 2026-03-08 (Build 9)
**Betrifft:** `symbolic_limit()`, `lhopital_applicable()`, `limit_comparison()`,
`partial_fraction_symbolic()`, `improper_integral_symbolic()`
**Beschreibung:**
`sympy.sympify()` wird verwendet, um Strings in SymPy-Ausdrücke zu konvertieren.
Intern verwendet sympify `eval()`, was bei unvertrautem Nutzer-Input zur
Ausführung beliebigen Python-Codes führen kann.
**Risikobewertung:** MEDIUM – Das Projekt ist eine lokale Python-Bibliothek
ohne Web-Interface. Die Funktionen werden nur von vertrauenswürdigen Entwicklern
aufgerufen. Kein direktes Produktionsrisiko.
**Empfehlung:** Falls das Projekt jemals als Web-Service genutzt wird:
`sympy.parsing.sympy_parser.parse_expr()` mit whitelist-basiertem Parser verwenden.
**Status:** Bekannt, akzeptiert (lokale Bibliothek, kein Web-Zugang)

---

## Sicherheitsanalyse (2026-03-08, Build 9)
Vollständige Überprüfung aller src/-Dateien:

| Risikokategorie | Status | Begründung |
|----------------|--------|------------|
| Code-Injection via sympify() | ⚠ MEDIUM | analysis.py – lokale Bibliothek, akzeptiert (s.o.) |
| Code-Injection (eval/exec direkt) | ✓ OK | Nicht direkt verwendet |
| Externe Prozesse (subprocess/os.system) | ✓ OK | Nicht verwendet |
| Unsichere Deserialisierung (pickle/yaml.load) | ✓ OK | Nicht verwendet |
| Netzwerkzugriffe | ✓ OK | Keine Netzwerkoperationen |
| Nutzereingaben (Webanwendung) | ✓ OK | Reine Python-Bibliothek, kein Web-Interface |
| Datenbankverbindungen | ✓ OK | Keine Datenbank |
| Kryptografisch unsichere Zufallszahlen (RSA) | ✓ OK | sympy.randprime() für Demo |
| Kryptografisch unsichere Zufallszahlen (Visualisierung) | ✓ OK | np.random.default_rng(42) nur für Sierpinski |
| Integer-Overflow | ✓ OK | Python hat unbegrenzte Integer |
| Division durch Null | ✓ OK | Alle relevanten Funktionen prüfen auf 0 |
| Path-Traversal (visualization.py save_path) | ✓ OK | Lokale Bibliothek, Nutzer kontrolliert Pfade selbst |
| ModularGroup Singularität (τ = -d/c) | ✓ OK | Liegt außerhalb des Fundamentalbereichs Im(τ)>0 |

## Neue Module Build 9 – Kurzprüfung

**modular_forms.py**: Keine externen Prozesse, keine Netzwerkzugriffe, keine Deserialisierung. Komplexe Arithmetik ohne Overflow-Risiken.

**p_adic.py**: Keine externen Prozesse. hensel_lift() nutzt modulare Arithmetik (kein Overflow). p_adic_exp/log haben Konvergenzprüfung.

**visualization.py**: matplotlib.use('Agg') korrekt gesetzt. save_path ist Nutzer-kontrolliert (akzeptiert für lokale Bibliothek).

## Hinweise
- Alle Funktionen validieren ihre Eingaben (ValueError bei ungültigen Parametern)
- RSA-Implementierung ist DEMO-only – für Produktion: PyCryptodome/cryptography verwenden
- Keine externen Abhängigkeiten mit bekannten CVEs (Stand 2026-03-08)
