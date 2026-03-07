# SECURITY_RISKS.md - specialist-maths

## Offene Sicherheitsrisiken
_(keine bekannten Risiken - geprüft am 2026-03-07)_

## Sicherheitsanalyse (2026-03-07)
Vollständige Überprüfung aller src/-Dateien:

| Risikokategorie | Status | Begründung |
|----------------|--------|------------|
| Code-Injection (eval/exec) | ✓ OK | Nicht verwendet |
| Externe Prozesse (subprocess/os.system) | ✓ OK | Nicht verwendet |
| Unsichere Deserialisierung (pickle/yaml.load) | ✓ OK | Nicht verwendet |
| Netzwerkzugriffe | ✓ OK | Keine Netzwerkoperationen |
| Nutzereingaben (Webanwendung) | ✓ OK | Reine Python-Bibliothek |
| Datenbankverbindungen | ✓ OK | Keine Datenbank |
| Kryptografisch unsichere Zufallszahlen | ✓ OK | random.Random nur für Monte-Carlo-Simulation (nicht sicherheitsrelevant) |
| Integer-Overflow | ✓ OK | Python hat unbegrenzte Integer |
| Division durch Null | ✓ OK | Wird in allen relevanten Funktionen abgefangen |

## Hinweise
- Alle Funktionen validieren ihre Eingaben (ValueError bei ungültigen Parametern)
- Keine externen Abhängigkeiten mit bekannten CVEs (Stand 2026-03-07)
- sympy, numpy, scipy sind aktuelle Versionen
