# SECURITY_FIXED.md - specialist-maths

## Behobene Sicherheitsprobleme

### sympify()-Code-Injection (Behoben 2026-03-09)
**Ursprünglich:** `sp.sympify()` in `analysis.py` verwendete intern `eval()`,
was bei unvertrauenswürdigen Nutzer-Strings zur Ausführung beliebigen Python-Codes
hätte führen können.

**Behebung:** Neue Hilfsfunktion `_safe_parse()` mit `parse_expr()` und
whitelist-basiertem `local_dict` (nur bekannte mathematische Funktionen erlaubt).
Fallback auf `sympify()` nur noch für SymPy-interne Ausdrücke (nicht Nutzereingaben).

**Betroffene Funktionen:**
- `symbolic_limit()`
- `lhopital_applicable()`
- `limit_comparison()`
- `partial_fraction_symbolic()`
- `improper_integral_symbolic()`

**Sicherheitsmechanismus:**
```python
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, \
    implicit_multiplication_application

def _safe_parse(expr_str, local_vars=None):
    transformations = standard_transformations + (implicit_multiplication_application,)
    safe_locals = {
        'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan,
        'exp': sp.exp, 'log': sp.log, ...
    }
    return parse_expr(expr_str, local_dict=safe_locals, transformations=transformations)
```

**Commit:** 2026-03-09, Build 10 (Vorbereitung)
