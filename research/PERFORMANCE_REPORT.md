# Performance-Report: specialist-maths

**Erstellt:** 2026-03-10
**Profiling-Tool:** Python `cProfile` + `time.perf_counter`
**Skript:** `/debugging/profile_performance.py`

---

## 1. Benchmark-Ergebnisse (Cold-Start vs. Warm-Start)

| Funktion                    | Aufrufe | Cold [ms] | Warm [ms] | Cold/Aufruf [µs] | Warm/Aufruf [µs] | Speedup |
|-----------------------------|---------|-----------|-----------|------------------|------------------|---------|
| `is_prime` (algebra_nt)     | 816     | 0.139     | 0.034     | 0.170            | 0.041            | 4.1×    |
| `prime_factorization`       | 356     | 0.145     | 0.015     | 0.408            | 0.041            | 9.9×    |
| `euler_phi`                 | 400     | 0.196     | 0.016     | 0.490            | 0.039            | 12.4×   |
| `bernoulli_number` (p_adic) | 70      | **1.060** | 0.003     | **15.148**       | 0.042            | **359.7×** |
| `p_adic_valuation`          | 846     | 0.181     | 0.047     | 0.214            | 0.056            | 3.8×    |
| `collatz_sequence`          | 178     | 0.173     | 0.007     | 0.973            | 0.038            | 25.5×   |
| `is_prime_fast` (proof)     | 786     | 0.139     | 0.030     | 0.176            | 0.039            | 4.6×    |

---

## 2. Cache-Statistiken (lru_cache)

| Funktion                | Hits   | Misses | Maxsize | Currsize |
|-------------------------|--------|--------|---------|----------|
| `is_prime`              | 1.134  | 498    | 10.000  | 498      |
| `prime_factorization`   | 510    | 202    | 1.000   | 202      |
| `euler_phi`             | 500    | 300    | 1.000   | 300      |
| `bernoulli_number`      | 120    | 20     | 200     | 20       |
| `p_adic_valuation`      | 1.095  | 597    | 5.000   | 597      |
| `collatz_sequence`      | 277    | 79     | 500     | 79       |
| `is_prime_fast`         | 1.074  | 498    | 10.000  | 498      |

---

## 3. Identifizierte Bottlenecks

### 3.1 Hauptbottleneck: `bernoulli_number` (p_adic.py)

**Größter Bottleneck** mit 15.15 µs pro Aufruf (Cold), 359.7× Cache-Speedup.

**Ursache laut cProfile:**
- 5.221 Funktionsaufrufe für nur 20 einzigartige Berechnungen
- Innere Funktion `binom()` (Binomialkoeffizienten) wird 2.660 Mal aufgerufen
- `min()` wird 2.470 Mal aufgerufen (innerhalb der Bernoulli-Schleife)
- Die Bernoulli-Formel `B_n = -Σ C(n+1,k) * B_k / (n+1)` ist O(n²) für B_n

**Bedeutung des Cache-Speedups von 359.7×:** Die Berechnung ist ohne Cache fast 360× langsamer. Das `lru_cache` kompensiert den Aufwand vollständig – sobald ein Wert berechnet ist, kostet der Abruf nur noch 0.042 µs.

### 3.2 Mittlerer Bottleneck: `euler_phi` (algebra_numbertheory.py)

**12.4× Cache-Speedup**, 0.49 µs pro Cold-Aufruf.

**Ursache:** `euler_phi(n)` ruft intern `prime_factorization(n)` auf, welche selbst gecacht ist. Aber der erste Aufruf muss die Faktorisierung berechnen (Probedivision bis √n).

### 3.3 Langsame erste Berechnung: `collatz_sequence`

**25.5× Cache-Speedup.** Collatz-Folgen können sehr lang sein (z.B. n=27 → 111 Schritte). Wiederholte Aufrufe werden vollständig durch den Cache abgedeckt.

### 3.4 Niedrige Bottleneck-Wirkung: `p_adic_valuation`

Trotz 846 Aufrufen ist der Speedup nur 3.8×, da die Basisfunktion bereits sehr effizient ist (einfache While-Schleife). Cache-Overhead pro Aufruf ist hier relativ hoch.

---

## 4. cProfile-Analyse der teuersten Funktionen (Cold-Start)

### `bernoulli_number` (20 unique Werte, 200 Aufrufe):
```
5221 function calls in 0.002 seconds
- bernoulli_number:  20 Aufrufe,  kumuliert 0.002s
- binom (intern):  2660 Aufrufe  ← HAUPTKOSTEN
- min (built-in):  2470 Aufrufe
```

### `euler_phi` (200 Aufrufe):
```
401 function calls in 0.000 seconds
- euler_phi: 200 Aufrufe, alle sehr schnell (gecachte Primfaktorisierung)
```

### `p_adic_valuation` (200 Aufrufe):
```
601 function calls in 0.000 seconds
- p_adic_valuation: 200 Aufrufe, abs() und isinstance() dominieren
- Sehr effizient: 0.214 µs/Aufruf Cold
```

---

## 5. Optimierungsempfehlungen

### Priorität HOCH

#### 5.1 `bernoulli_number`: Verwende Zolotarev-Formel oder mpmath
**Problem:** Aktuelle O(n²) Implementierung via Bernoulli-Rekursion.
**Lösung A:** `mpmath.bernoulli(n)` verwenden (hochoptimierte C-Bibliothek).
**Lösung B:** Vorberechnete Tabelle für B_0 bis B_50 hardcoden (häufig benötigte Werte).
**Speedup-Schätzung:** 10–50× für große n.

#### 5.2 `binom`-Hilfsfunktion in `bernoulli_number`: Caching hinzufügen
**Problem:** `binom()` wird innerhalb von `bernoulli_number` wiederholt mit denselben Argumenten aufgerufen (2.660 Aufrufe für 20 unique Bernoulli-Zahlen).
**Lösung:** `@functools.lru_cache` auf die interne `binom()`-Funktion anwenden.
**Speedup-Schätzung:** 2–5× für Cold-Starts.

### Priorität MITTEL

#### 5.3 `eisenstein_series`: Symmetrie ausnutzen (bereits implementiert als `eisenstein_series_fast`)
**Status:** Implementiert in Build 13.
**Erwarteter Speedup:** ~4× für gerades k ≥ 4.

#### 5.4 `goldbach_verification_range`: Multiprocessing (bereits implementiert)
**Status:** Implementiert in Build 13 mit `multiprocessing.Pool`.
**Erwarteter Speedup:** 2–4× auf Mehrkern-Systemen.

#### 5.5 `prime_factorization`: Pollard-Rho-Algorithmus für große n
**Problem:** Aktuelle Probedivision ist O(√n).
**Lösung:** Pollard-Rho für n > 10^6 (O(n^{1/4})).
**Impact:** Nur für sehr große Zahlen relevant; für n < 10^6 ist Probedivision optimal.

### Priorität NIEDRIG

#### 5.6 `p_adic_valuation`: isinstance()-Aufruf entfernen
**Problem:** `isinstance()` wird bei jedem der 846 Aufrufe ausgeführt.
**Lösung:** Typprüfung nur einmalig oder als Assertion-Mode implementieren.
**Speedup:** < 10%, Aufwand minimal.

#### 5.7 Vektorialisierung der Collatz-Berechnung mit NumPy
**Problem:** Reine Python-Schleife für die Collatz-Folge.
**Lösung:** NumPy-vektorisierte Version für Batch-Berechnungen.
**Limitierung:** Collatz-Folgen haben unterschiedliche Längen → NumPy-Padding nötig.

---

## 6. Architekturelle Empfehlungen

### 6.1 Globaler Primzahl-Cache (Sieve)
Anstatt `lru_cache` pro Funktion: **Einmaliger Sieb-Cache** bei Modulinitialisierung.
Alle primzahlabhängigen Funktionen greifen auf denselben Sieb zurück.
Speicherkosten: ~1 MB für alle Primes bis 10^7.

### 6.2 Lazy Initialization für teure Berechnungen
`bernoulli_number`, `euler_phi` etc. könnten als **class-level cached properties** implementiert werden, um Neuberechnungen bei Modul-Reloads zu vermeiden.

### 6.3 C-Extension für kritische Pfade
Falls Performance-kritisch: `is_prime()` und `prime_factorization()` als Cython- oder CFFI-Extension (10–100× Speedup möglich).

---

## 7. Zusammenfassung

**Stärken:**
- `lru_cache` ist gut eingesetzt und liefert bis zu 360× Speedup
- Cache-Hit-Raten sind hoch (z.B. 1.134 Hits bei `is_prime`)
- Alle Funktionen terminieren in Microsekunden (Warm)

**Schwächen:**
- `bernoulli_number` ist der einzige signifikante Bottleneck (15 µs/Aufruf Cold)
- Kein vektorisiertes Batch-Interface für Primzahl-Berechnungen
- `binom()` innerhalb `bernoulli_number` nicht gecacht

**Gesamteinschätzung:** Das Modul ist für wissenschaftliche Berechnungen gut optimiert. Kritischer Einsatz (z.B. n > 10^6, tausende Bernoulli-Zahlen) würde `mpmath` oder C-Extensions erfordern.
