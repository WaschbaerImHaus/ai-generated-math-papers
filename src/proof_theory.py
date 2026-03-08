"""
Beweistheorie-Modul für mathematische Beweise und Vermutungsanalyse.

Dieses Modul implementiert Werkzeuge zur algorithmischen Untersuchung
offener mathematischer Vermutungen. Es dient als Grundlage für das
Fernziel, bedeutende mathematische Vermutungen zu beweisen oder zu widerlegen.

@author: Kurt Ingwer
@version: 1.0
@since: 2026-03-08
@lastModified: 2026-03-08
"""

import math
import itertools
from typing import Callable, Generator, Optional
import sympy
import numpy as np


# ===========================================================================
# COLLATZ-VERMUTUNG (3n+1-Problem)
# ===========================================================================

def collatz_sequence(n: int) -> list[int]:
    """
    Berechnet die Collatz-Folge für eine gegebene Startzahl.

    Die Collatz-Vermutung besagt, dass diese Folge für jede positive
    ganze Zahl n nach endlich vielen Schritten die Zahl 1 erreicht.

    Regel:
        f(n) = n/2      falls n gerade
        f(n) = 3n + 1   falls n ungerade

    @param n: Startzahl (positive ganze Zahl)
    @return: Liste aller Glieder der Folge bis 1 (inkl. n und 1)
    @raises ValueError: Wenn n ≤ 0
    @lastModified: 2026-03-08
    """
    if n <= 0:
        raise ValueError(f"n muss positiv sein, erhalten: {n}")

    sequence = [n]
    current = n
    while current != 1:
        if current % 2 == 0:
            current = current // 2      # gerade: halbieren
        else:
            current = 3 * current + 1  # ungerade: 3n+1
        sequence.append(current)

    return sequence


def collatz_stopping_time(n: int) -> int:
    """
    Gibt die Anzahl der Schritte zurück, bis die Collatz-Folge 1 erreicht.

    Diese "Stoppzeit" ist ein Maß für die Komplexität der Zahl n.
    Zahlen mit langen Stoppzeiten sind oft Kandidaten für Musteranalysen.

    @param n: Startzahl (positive ganze Zahl)
    @return: Anzahl der Schritte bis 1
    @lastModified: 2026-03-08
    """
    steps = 0
    current = n
    while current != 1:
        if current % 2 == 0:
            current = current // 2
        else:
            current = 3 * current + 1
        steps += 1
    return steps


def collatz_max_value(n: int) -> int:
    """
    Gibt den Maximalwert zurück, den die Collatz-Folge bei Startzahl n erreicht.

    @param n: Startzahl
    @return: Maximaler Wert in der Folge
    @lastModified: 2026-03-08
    """
    return max(collatz_sequence(n))


def collatz_verify_range(limit: int) -> dict:
    """
    Verifiziert die Collatz-Vermutung für alle Zahlen von 1 bis limit.

    Dies ist eine empirische Überprüfung – kein formaler Beweis.
    Bekannt: Verifiziert bis ca. 2^68 ≈ 2.95 × 10^20.

    @param limit: Obere Grenze der Verifikation (inklusiv)
    @return: Dictionary mit Statistiken:
             - verified: bool (alle konvergiert)
             - counterexample: Optional[int] (erste Gegenbeispiel oder None)
             - max_stopping_time: (zahl, schritte) mit längster Stoppzeit
             - max_value_seen: (zahl, maxwert) mit größtem Zwischenwert
    @lastModified: 2026-03-08
    """
    max_steps = (0, 0)          # (startzahl, schritte)
    max_value = (0, 0)          # (startzahl, maxwert)
    counterexample = None

    for n in range(1, limit + 1):
        try:
            seq = collatz_sequence(n)
            steps = len(seq) - 1
            peak = max(seq)

            if steps > max_steps[1]:
                max_steps = (n, steps)
            if peak > max_value[1]:
                max_value = (n, peak)
        except RecursionError:
            counterexample = n
            break

    return {
        "verified": counterexample is None,
        "counterexample": counterexample,
        "max_stopping_time": max_steps,
        "max_value_seen": max_value,
        "range_checked": limit
    }


# ===========================================================================
# GOLDBACH-VERMUTUNG
# ===========================================================================

def is_prime_fast(n: int) -> bool:
    """
    Schneller deterministischer Primzahltest für mittlere Zahlen.

    Verwendet die 6k±1-Optimierung: Alle Primzahlen > 3 haben die Form
    6k+1 oder 6k-1, weil alle anderen durch 2 oder 3 teilbar sind.

    @param n: Zu prüfende Zahl
    @return: True wenn prim, False sonst
    @lastModified: 2026-03-08
    """
    if n < 2:
        return False
    if n < 4:
        return True             # 2, 3 sind prim
    if n % 2 == 0 or n % 3 == 0:
        return False            # Teiler 2 oder 3

    # Teste nur Zahlen der Form 6k±1 bis √n
    k = 5
    while k * k <= n:
        if n % k == 0 or n % (k + 2) == 0:
            return False
        k += 6

    return True


def goldbach_decomposition(n: int) -> Optional[tuple[int, int]]:
    """
    Findet eine Goldbach-Zerlegung der geraden Zahl n in zwei Primzahlen.

    Die Goldbach-Vermutung (1742) besagt, dass jede gerade Zahl > 2
    als Summe zweier Primzahlen darstellbar ist.

    Strategie: Für p von 2 bis n/2 prüfen ob p prim und n-p prim.

    @param n: Gerade ganze Zahl > 2
    @return: Tupel (p, q) mit p + q = n und p, q prim, oder None
    @raises ValueError: Wenn n ungerade oder n ≤ 2
    @lastModified: 2026-03-08
    """
    if n <= 2 or n % 2 != 0:
        raise ValueError(f"n muss gerade und > 2 sein, erhalten: {n}")

    # Suche die "schönste" Zerlegung nahe n/2
    for p in range(2, n // 2 + 1):
        q = n - p
        if is_prime_fast(p) and is_prime_fast(q):
            return (p, q)

    return None  # Sollte nie passieren, wenn Goldbach gilt


def goldbach_all_decompositions(n: int) -> list[tuple[int, int]]:
    """
    Gibt ALLE Goldbach-Zerlegungen von n zurück.

    @param n: Gerade ganze Zahl > 2
    @return: Liste aller Paare (p, q) mit p ≤ q und p + q = n
    @lastModified: 2026-03-08
    """
    if n <= 2 or n % 2 != 0:
        raise ValueError(f"n muss gerade und > 2 sein, erhalten: {n}")

    decompositions = []
    for p in range(2, n // 2 + 1):
        q = n - p
        if is_prime_fast(p) and is_prime_fast(q):
            decompositions.append((p, q))

    return decompositions


def goldbach_verify_range(limit: int) -> dict:
    """
    Verifiziert die Goldbach-Vermutung für alle geraden Zahlen bis limit.

    @param limit: Obere Grenze (wird auf nächste gerade Zahl abgerundet)
    @return: Dictionary mit Statistiken
    @lastModified: 2026-03-08
    """
    verified = True
    counterexample = None
    min_decompositions = (4, 1)  # (zahl, anzahl) – 4 hat nur eine Zerlegung

    for n in range(4, limit + 1, 2):
        decomps = goldbach_all_decompositions(n)
        if not decomps:
            verified = False
            counterexample = n
            break
        if len(decomps) < min_decompositions[1] or min_decompositions[0] == 4:
            min_decompositions = (n, len(decomps))

    return {
        "verified": verified,
        "counterexample": counterexample,
        "min_decompositions": min_decompositions,
        "range_checked": limit
    }


# ===========================================================================
# ZWILLINGSPRIMZAHL-VERMUTUNG
# ===========================================================================

def find_twin_primes(limit: int) -> list[tuple[int, int]]:
    """
    Findet alle Zwillingsprimzahlpaare (p, p+2) bis zur Schranke limit.

    Zwillingsprimzahlpaare: (3,5), (5,7), (11,13), (17,19), (29,31), ...
    Die Vermutung: Es gibt unendlich viele solche Paare.

    @param limit: Obere Schranke
    @return: Liste aller Zwillingsprimzahlpaare (p, p+2) mit p+2 ≤ limit
    @lastModified: 2026-03-08
    """
    twins = []
    for p in range(3, limit - 1):
        if is_prime_fast(p) and is_prime_fast(p + 2):
            twins.append((p, p + 2))
    return twins


def twin_prime_count(limit: int) -> int:
    """
    Zählt Zwillingsprimzahlpaare bis zur Schranke limit.

    Die Hardy-Littlewood-Vermutung gibt eine asymptotische Schätzung:
    π₂(x) ~ 2·C₂ · x / (ln x)²
    mit der Zwillingsprimzahlkonstante C₂ ≈ 0.6601618...

    @param limit: Obere Schranke
    @return: Anzahl der Zwillingsprimzahlpaare
    @lastModified: 2026-03-08
    """
    return len(find_twin_primes(limit))


def twin_prime_constant() -> float:
    """
    Approximiert die Zwillingsprimzahlkonstante C₂ ≈ 0.6601618...

    C₂ = Π_{p prim, p>2} p(p-2)/(p-1)²

    Diese Konstante erscheint in der Hardy-Littlewood-Vermutung für
    die Anzahl der Zwillingsprimzahlpaare bis x.

    @return: Näherungswert von C₂
    @lastModified: 2026-03-08
    """
    # Produkt über die ersten 1000 ungeraden Primzahlen
    product = 1.0
    count = 0
    n = 3
    while count < 1000:
        if is_prime_fast(n):
            product *= n * (n - 2) / ((n - 1) ** 2)
            count += 1
        n += 2
    return product


# ===========================================================================
# RIEMANN-ZETA-FUNKTION
# ===========================================================================

def riemann_zeta_partial(s: complex, terms: int = 10000) -> complex:
    """
    Berechnet die Riemann-Zeta-Funktion via Partialsumme für Re(s) > 1.

    ζ(s) = Σ(n=1 bis ∞) 1/n^s

    WARNUNG: Konvergiert nur für Re(s) > 1. Für analytische Fortsetzung
    wird die Euler-Maclaurin-Formel oder die Euler-Alternating-Series benötigt.

    @param s: Komplexe Zahl s mit Re(s) > 1
    @param terms: Anzahl der Summanden (mehr = genauer, aber langsamer)
    @return: Näherungswert von ζ(s)
    @raises ValueError: Wenn Re(s) ≤ 1
    @lastModified: 2026-03-08
    """
    if s.real <= 1:
        raise ValueError(
            f"Partialsumme konvergiert nur für Re(s) > 1, erhalten: Re(s) = {s.real}"
        )

    total = complex(0.0)
    for n in range(1, terms + 1):
        total += 1.0 / (n ** s)
    return total


def riemann_zeta_euler_maclaurin(s: complex, terms: int = 100) -> complex:
    """
    Berechnet ζ(s) via Euler-Maclaurin-Formel – konvergiert schneller.

    Verwendet die Beschleunigung durch Euler-Maclaurin:
    ζ(s) ≈ Σ(n=1 bis N) 1/n^s + N^(1-s)/(s-1) + N^(-s)/2 + ...

    Funktioniert für Re(s) > 0 (außer s=1).

    @param s: Komplexe Zahl s (Re(s) > 0, s ≠ 1)
    @param terms: Anzahl der direkten Summanden vor Euler-Maclaurin-Korrektur
    @return: Näherungswert von ζ(s)
    @lastModified: 2026-03-08
    """
    if abs(s - 1) < 1e-10:
        raise ValueError("ζ(s) hat einen Pol bei s = 1")

    # Partialsumme bis N
    N = terms
    total = complex(0.0)
    for n in range(1, N + 1):
        total += 1.0 / (n ** s)

    # Euler-Maclaurin-Hauptkorrektur: Restintegral ∫_N^∞ x^{-s} dx = N^{1-s}/(s-1)
    total += (N ** (1 - s)) / (s - 1)

    # Erster Euler-Maclaurin-Term: f(N)/2
    total += 0.5 / (N ** s)

    # Bernoulli-Korrektur (B_2/2! · f'(N), B_4/4! · f'''(N), ...)
    # B_2 = 1/6, Ableitung von n^{-s} ist -s·n^{-s-1}
    total += (1.0 / 12.0) * (-s) / (N ** (s + 1))
    # B_4 = -1/30
    total += (-1.0 / 720.0) * (-s) * (-s - 1) * (-s - 2) / (N ** (s + 3))

    return total


def check_riemann_hypothesis_numerically(
    imaginary_range: tuple[float, float],
    resolution: int = 100
) -> dict:
    """
    Numerische Überprüfung der Riemann-Hypothese im angegebenen Bereich.

    Sucht Nullstellen von ζ(1/2 + it) für t im angegebenen Bereich.
    Die RH behauptet, dass ALLE nicht-trivialen Nullstellen auf Re(s) = 1/2 liegen.

    Bekannte Nullstellen (Imaginärteile):
    t₁ ≈ 14.1347, t₂ ≈ 21.0220, t₃ ≈ 25.0109, t₄ ≈ 30.4249, ...

    @param imaginary_range: (t_min, t_max) – Bereich des Imaginärteils
    @param resolution: Anzahl der Testpunkte
    @return: Dictionary mit gefundenen Nullstellen-Näherungen
    @lastModified: 2026-03-08
    """
    t_min, t_max = imaginary_range
    t_values = np.linspace(t_min, t_max, resolution)

    # Zeta-Werte auf der kritischen Geraden berechnen
    zeta_values = []
    for t in t_values:
        s = complex(0.5, t)
        try:
            z = riemann_zeta_euler_maclaurin(s, terms=200)
            zeta_values.append(abs(z))
        except Exception:
            zeta_values.append(float('inf'))

    # Lokale Minima finden (Nullstellenkandidaten)
    zero_candidates = []
    for i in range(1, len(zeta_values) - 1):
        if (zeta_values[i] < zeta_values[i-1] and
                zeta_values[i] < zeta_values[i+1] and
                zeta_values[i] < 0.5):  # Nahe Null
            zero_candidates.append({
                "t": t_values[i],
                "abs_zeta": zeta_values[i],
                "on_critical_line": True  # Per Konstruktion
            })

    return {
        "range": imaginary_range,
        "resolution": resolution,
        "zero_candidates": zero_candidates,
        "known_zeros_in_range": [
            t for t in [14.1347, 21.0220, 25.0109, 30.4249, 32.9351,
                        37.5862, 40.9187, 43.3271, 48.0052, 49.7738]
            if t_min <= t <= t_max
        ]
    }


# ===========================================================================
# FORMALE BEWEISTECHNIKEN
# ===========================================================================

class ProofByInduction:
    """
    Hilfsklasse für den Aufbau von Beweisen durch vollständige Induktion.

    Schema:
        1. Basisfall: Zeige P(n₀) für den kleinsten Fall n₀
        2. Induktionsschritt: Zeige P(n) → P(n+1)
        3. Schluss: P(n) gilt für alle n ≥ n₀

    @author: Kurt Ingwer
    @since: 2026-03-08
    @lastModified: 2026-03-08
    """

    @staticmethod
    def verify_base_case(predicate: Callable[[int], bool], base: int) -> bool:
        """
        Überprüft den Basisfall einer Induktion.

        @param predicate: Die zu beweisende Aussage P(n) als Python-Funktion
        @param base: Der Basisfall n₀
        @return: True wenn P(base) gilt
        @lastModified: 2026-03-08
        """
        return predicate(base)

    @staticmethod
    def empirical_verify(
        predicate: Callable[[int], bool],
        start: int,
        end: int
    ) -> dict:
        """
        Empirische Überprüfung einer Aussage für einen Bereich.

        WICHTIG: Dies ist KEIN Beweis! Empirische Verifikation kann nur
        Gegenbeispiele finden, aber keine Aussage beweisen.

        @param predicate: Die Aussage P(n)
        @param start: Anfang des Bereichs
        @param end: Ende des Bereichs (inklusiv)
        @return: Verifikationsergebnis mit eventuellem Gegenbeispiel
        @lastModified: 2026-03-08
        """
        for n in range(start, end + 1):
            if not predicate(n):
                return {
                    "verified": False,
                    "counterexample": n,
                    "range": (start, end),
                    "warning": "Gegenbeispiel gefunden – Aussage ist FALSCH"
                }
        return {
            "verified": True,
            "counterexample": None,
            "range": (start, end),
            "warning": "Empirisch bestätigt – KEIN formaler Beweis!"
        }


def sieve_of_eratosthenes(limit: int) -> list[int]:
    """
    Sieb des Eratosthenes – findet alle Primzahlen bis limit.

    Algorithmus (ca. 240 v. Chr.):
        1. Erstelle Boole-Array von 2 bis limit, alle True
        2. Für jede Primzahl p: markiere alle Vielfachen p², p²+p, ... als False
        3. Übrige True-Einträge sind Primzahlen

    Zeitkomplexität: O(n log log n)
    Speicherkomplexität: O(n)

    @param limit: Obere Grenze (inklusiv)
    @return: Liste aller Primzahlen ≤ limit
    @lastModified: 2026-03-08
    """
    if limit < 2:
        return []

    # Boole-Array: is_prime[i] = True bedeutet i ist prim
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False

    # Sieben: alle Vielfachen von p markieren, beginnend bei p²
    p = 2
    while p * p <= limit:
        if is_prime[p]:
            # Alle Vielfachen von p ab p² als nicht-prim markieren
            for multiple in range(p * p, limit + 1, p):
                is_prime[multiple] = False
        p += 1

    return [n for n in range(2, limit + 1) if is_prime[n]]


def miller_rabin_primality_test(n: int, rounds: int = 20) -> bool:
    """
    Miller-Rabin-Primzahltest – probabilistisch, für sehr große Zahlen.

    Für rounds ≥ 20 ist die Fehlerwahrscheinlichkeit < 4^(-20) ≈ 10^(-12).
    Mit den festen Basen {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}
    ist der Test deterministisch für alle n < 3.317 × 10^24.

    Grundlage: Fermatscher Kleiner Satz – wenn p prim, dann a^(p-1) ≡ 1 (mod p).

    @param n: Zu prüfende Zahl (n ≥ 2)
    @param rounds: Anzahl der Zeugen-Tests
    @return: True wenn wahrscheinlich prim, False wenn definitiv zusammengesetzt
    @lastModified: 2026-03-08
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    # Schreibe n-1 = 2^r · d mit d ungerade
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # Deterministische Basen für n < 3.317 × 10^24
    deterministic_bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    test_bases = [b for b in deterministic_bases if b < n]

    # Wenn n sehr klein, fallback auf trial division
    if not test_bases:
        return is_prime_fast(n)

    for a in test_bases:
        # Berechne a^d mod n
        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue  # Kein Zeuge für Zusammengesetztheit

        # Quadriere r-1 mal
        composite = True
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                composite = False
                break

        if composite:
            return False  # Definitiv zusammengesetzt

    return True  # Wahrscheinlich prim (für obige Basen: deterministisch)


# ===========================================================================
# ZAHLENTHEORETISCHE HILFSFUNKTIONEN (für Beweise)
# ===========================================================================

def legendre_symbol(a: int, p: int) -> int:
    """
    Berechnet das Legendre-Symbol (a/p) für ungerade Primzahl p.

    (a/p) = 0  falls p | a
    (a/p) = 1  falls a quadratischer Rest mod p
    (a/p) = -1 falls a quadratischer Nichtrest mod p

    Berechnung via Euler-Kriterium: (a/p) ≡ a^((p-1)/2) (mod p)

    @param a: Ganzzahl
    @param p: Ungerade Primzahl
    @return: Legendre-Symbol ∈ {-1, 0, 1}
    @raises ValueError: Wenn p keine ungerade Primzahl ist
    @lastModified: 2026-03-08
    """
    if not miller_rabin_primality_test(p) or p == 2:
        raise ValueError(f"p muss eine ungerade Primzahl sein, erhalten: {p}")

    a = a % p
    if a == 0:
        return 0

    # Euler-Kriterium
    result = pow(a, (p - 1) // 2, p)
    return -1 if result == p - 1 else result


def jacobi_symbol(a: int, n: int) -> int:
    """
    Berechnet das Jacobi-Symbol (a/n) für ungerades n > 0.

    Verallgemeinerung des Legendre-Symbols auf zusammengesetzte Zahlen.
    Verwendet die quadratischen Reziprozitätsgesetze für effiziente Berechnung.

    WICHTIG: (a/n) = 1 bedeutet NICHT, dass a quadratischer Rest mod n ist!

    @param a: Ganzzahl
    @param n: Ungerade positive ganze Zahl > 1
    @return: Jacobi-Symbol ∈ {-1, 0, 1}
    @lastModified: 2026-03-08
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError(f"n muss ungerade und positiv sein, erhalten: {n}")

    a = a % n
    result = 1

    while a != 0:
        # Teile a durch 2, solange möglich
        while a % 2 == 0:
            a //= 2
            # Quadratisches Reziprositätsgesetz für 2:
            # (2/n) = (-1)^((n²-1)/8)
            if n % 8 in (3, 5):
                result = -result

        # Reziprozität: (a/n)(n/a) = (-1)^((a-1)(n-1)/4)
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result

        a = a % n

    return result if n == 1 else 0


def chinese_remainder_theorem(remainders: list[int], moduli: list[int]) -> int:
    """
    Chinesischer Restsatz (CRT): Löst simultane Kongruenzen.

    Gegeben: x ≡ r₁ (mod m₁), x ≡ r₂ (mod m₂), ..., x ≡ rₖ (mod mₖ)
    Gesucht: Kleinste nichtnegative Lösung x, falls die mᵢ paarweise koprim sind.

    Konstruktion via Bézout:
        M = m₁ · m₂ · ... · mₖ
        Mᵢ = M / mᵢ
        yᵢ = Mᵢ⁻¹ mod mᵢ (Modulares Inverses)
        x = Σ rᵢ · Mᵢ · yᵢ (mod M)

    @param remainders: Liste der Reste [r₁, r₂, ..., rₖ]
    @param moduli: Liste der Moduln [m₁, m₂, ..., mₖ] (paarweise koprim)
    @return: Kleinste nichtnegative Lösung x
    @raises ValueError: Wenn Moduln nicht paarweise koprim
    @lastModified: 2026-03-08
    """
    if len(remainders) != len(moduli):
        raise ValueError("Anzahl der Reste muss mit Anzahl der Moduln übereinstimmen")

    # Paarweise Koprimheit prüfen
    for i in range(len(moduli)):
        for j in range(i + 1, len(moduli)):
            if math.gcd(moduli[i], moduli[j]) != 1:
                raise ValueError(
                    f"Moduln müssen paarweise koprim sein: "
                    f"gcd({moduli[i]}, {moduli[j]}) = {math.gcd(moduli[i], moduli[j])}"
                )

    M = 1
    for m in moduli:
        M *= m             # Gesamtmodul M = m₁ · m₂ · ... · mₖ

    x = 0
    for r, m in zip(remainders, moduli):
        Mi = M // m        # Mᵢ = M / mᵢ
        # Modulares Inverses von Mᵢ mod mᵢ (via Kleinem Fermatschen Satz wenn mᵢ prim)
        yi = pow(Mi, -1, m)
        x += r * Mi * yi

    return x % M


# ===========================================================================
# BEWEISSTATISTIKEN UND -DOKUMENTATION
# ===========================================================================

# ===========================================================================
# KREISMETHODE (HARDY-LITTLEWOOD) FÜR GOLDBACH-VERMUTUNG
# ===========================================================================

def hardy_littlewood_singular_series(n: int) -> float:
    """
    Berechnet die Singuläre Reihe S(n) aus der Hardy-Littlewood-Kreismethode.

    Die singuläre Reihe entsteht als Eulerprodukt und codiert die lokalen
    (p-adischen) Lösungsdichten der Goldbach-Gleichung p₁ + p₂ = n.

    Für gerades n gilt:
        S(n) = Π_{p | n, p > 2} [p(p-2)/(p-1)²]⁻¹ · p(p-2)/(p-1)²
             · Π_{p ∤ n, p > 2} p(p-2)/(p-1)²  /  Π_{p > 2} p(p-2)/(p-1)²

    Vereinfacht (Hardy-Littlewood, 1923):
        S(n) = Π_{p | n, p prim > 2} (p-1)/(p-2)
             · Π_{p ∤ n, p prim > 2} 1   (Beitrag = 1 für p ∤ n)

    D.h. der Faktor ≠ 1 kommt nur von Primteilern p > 2 von n.

    @param n: Gerade positive ganze Zahl ≥ 4
    @return: Wert der Singulären Reihe S(n) ≥ 1
    @raises ValueError: Wenn n ungerade oder n < 4
    @lastModified: 2026-03-08
    """
    if n < 4 or n % 2 != 0:
        raise ValueError(f"n muss gerade und ≥ 4 sein, erhalten: {n}")

    # Ermittle alle Primteiler von n, die größer als 2 sind
    # (p=2 trägt zum Eulerprodukt separat bei – für gerades n immer 2|n)
    singular_series = 1.0
    temp = n
    # Entferne Faktoren 2 (die nicht zur Singulären Reihe beitragen für gerades n)
    while temp % 2 == 0:
        temp //= 2

    # Sammle ungerade Primteiler von n
    odd_prime_factors = set()
    d = 3
    t = temp
    while d * d <= t:
        if t % d == 0:
            odd_prime_factors.add(d)
            while t % d == 0:
                t //= d
        d += 2
    if t > 1:
        odd_prime_factors.add(t)

    # Für jeden ungeraden Primteiler p von n: Faktor (p-1)/(p-2)
    # Dieser Faktor ist > 1 und erhöht die Anzahl der Zerlegungen
    # (mehr Möglichkeiten, weil p | n erlaubt p als direkten Summanden)
    for p in odd_prime_factors:
        if p > 2:
            singular_series *= (p - 1) / (p - 2)

    return singular_series


def goldbach_circle_method_estimate(n: int) -> float:
    """
    Schätzt die Anzahl der Goldbach-Zerlegungen via Hardy-Littlewood-Kreismethode.

    Die Kreismethode (Hardy & Littlewood, 1923) liefert eine asymptotische Formel:
        r(n) ≈ 2 · C₂ · n / (ln n)² · S(n)

    wobei:
        C₂ = Π_{p prim, p > 2} p(p-2)/(p-1)² ≈ 0.6601618...
             die Zwillingsprimzahlkonstante ist
        S(n) = singuläre Reihe, die lokale Arithmetik codiert

    Interpretation: Für großes gerades n gibt es asymptotisch
    ca. r(n) Zerlegungen n = p + q mit p, q prim.

    @param n: Gerade positive ganze Zahl ≥ 4
    @return: Geschätzte Anzahl der Goldbach-Zerlegungen (ungeordnet, p ≤ q)
    @raises ValueError: Wenn n ungerade oder n < 4
    @lastModified: 2026-03-08
    """
    if n < 4 or n % 2 != 0:
        raise ValueError(f"n muss gerade und ≥ 4 sein, erhalten: {n}")

    # Zwillingsprimzahlkonstante C₂ (bekannter Wert)
    twin_prime_constant_C2 = 0.6601618158468695

    # Singuläre Reihe berechnen
    S_n = hardy_littlewood_singular_series(n)

    # Hardy-Littlewood-Asymptotik: r(n) ≈ 2 · C₂ · n / (ln n)² · S(n)
    # Der Faktor 2 kommt davon, dass die Goldbach-Gleichung p+q=n symmetrisch ist
    ln_n = math.log(n)
    estimate = 2.0 * twin_prime_constant_C2 * n / (ln_n ** 2) * S_n

    # Die Schätzung gilt für UNGEORDNETE Paare (p ≤ q), daher halbieren
    return estimate / 2.0


def goldbach_circle_method_accuracy(n_max: int) -> list[dict]:
    """
    Vergleicht die Hardy-Littlewood-Schätzung mit tatsächlichen Zerlegungszahlen.

    Für jede gerade Zahl n von 4 bis n_max wird berechnet:
        - Die tatsächliche Anzahl der Goldbach-Zerlegungen (exakt)
        - Die Schätzung via Kreismethode
        - Das Verhältnis actual / estimate (sollte → 1 für großes n)

    @param n_max: Obere Grenze (wird auf nächste gerade Zahl abgerundet)
    @return: Liste von Dictionaries mit Feldern:
             'n': die gerade Zahl
             'actual': tatsächliche Zerlegungsanzahl
             'estimate': Schätzung via Kreismethode
             'ratio': Verhältnis actual/estimate (Güte der Approximation)
    @lastModified: 2026-03-08
    """
    results = []

    for n in range(4, n_max + 1, 2):
        # Tatsächliche Anzahl der Zerlegungen (exakt gezählt, p ≤ q)
        actual_count = len(goldbach_all_decompositions(n))

        # Schätzung via Kreismethode
        estimate = goldbach_circle_method_estimate(n)

        # Verhältnis (Güte der Approximation, → 1 für n → ∞)
        ratio = actual_count / estimate if estimate > 0 else float('inf')

        results.append({
            'n': n,
            'actual': actual_count,
            'estimate': estimate,
            'ratio': ratio
        })

    return results


def conjecture_status_report() -> dict:
    """
    Gibt eine Übersicht über alle bekannten offenen Vermutungen zurück.

    @return: Dictionary mit Status aller bekannten Vermutungen
    @lastModified: 2026-03-08
    """
    return {
        "millennium_problems": {
            "Riemann-Hypothese": {
                "year": 1859,
                "status": "offen",
                "verified_up_to": "10^13 Nullstellen",
                "prize_usd": 1_000_000
            },
            "P_vs_NP": {
                "year": 1971,
                "status": "offen",
                "prize_usd": 1_000_000
            },
            "Hodge-Vermutung": {
                "year": 1950,
                "status": "offen",
                "prize_usd": 1_000_000
            },
            "Yang-Mills": {
                "year": 2000,
                "status": "offen",
                "prize_usd": 1_000_000
            },
            "Navier-Stokes": {
                "year": 2000,
                "status": "offen",
                "prize_usd": 1_000_000
            },
            "BSD-Vermutung": {
                "year": 1965,
                "status": "offen",
                "prize_usd": 1_000_000
            },
            "Poincaré-Vermutung": {
                "year": 1904,
                "status": "bewiesen (Perelman 2003)",
                "prize_usd": 1_000_000
            }
        },
        "other_conjectures": {
            "Goldbach-Vermutung": {
                "year": 1742,
                "status": "offen",
                "verified_up_to": "4 × 10^18"
            },
            "Collatz-Vermutung": {
                "year": 1937,
                "status": "offen",
                "verified_up_to": "≈ 2^68"
            },
            "Zwillingsprimzahl": {
                "year": "Antike",
                "status": "offen",
                "partial": "Zhang 2013: Abstände < 246 unendlich oft"
            },
            "ABC-Vermutung": {
                "year": 1985,
                "status": "umstrittener Beweis (Mochizuki 2012)"
            }
        }
    }
