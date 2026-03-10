"""
@file algebra_diophantine.py
@brief Diophantische Gleichungen, quadratische Reste, Tonelli-Shanks, Cipolla.
@description
    Enthält Algorithmen für ganzzahlige Gleichungen und modulare Quadratwurzeln:

    Diophantische Gleichungen:
    - solve_linear_diophantine()          – ax + by = c via Bezout
    - solve_quadratic_diophantine_pell()  – Pell-Gleichung x² - Dy² = 1
    - solve_pythagorean_triples()         – Pythagoräische Tripel (Euklid's Formel)
    - solve_diophantine_two_squares()     – n = a² + b² (Fermat's Zwei-Quadrate-Satz)
    - markov_numbers()                    – Markov-Zahlen (Markov-Gleichung)

    Quadratische Reste und Reziprozität:
    - is_quadratic_residue()   – Euler-Kriterium
    - quadratic_residues()     – Alle QR modulo p
    - quadratic_reciprocity()  – Gauß'sches Reziprozitätsgesetz

    Quadratwurzeln modulo Primzahlen:
    - tonelli_shanks()         – √n mod p (allgemeiner Algorithmus)
    - cipolla_algorithm()      – √n mod p (Cipolla 1903, via GF(p²))

    Importiert gcd, extended_gcd, mod_inverse, is_prime aus den Kernmodulen.
    Ausgelagert aus algebra.py für bessere Modularität.

@author Kurt Ingwer
@lastModified 2026-03-10
"""

import math

# Abhängigkeiten aus den Kernmodulen importieren
from algebra_core import gcd, extended_gcd, mod_inverse
from algebra_numbertheory import is_prime


# ===========================================================================
# DIOPHANTISCHE GLEICHUNGEN
# ===========================================================================

def solve_linear_diophantine(a: int, b: int, c: int) -> tuple[int, int, int] | None:
    """
    @brief Löst die lineare Diophantische Gleichung ax + by = c in ganzen Zahlen.
    @description
        Eine lineare Diophantische Gleichung hat ganzzahlige Lösungen genau dann,
        wenn gcd(a, b) | c (d.h. gcd(a,b) teilt c).

        Algorithmus:
        1. g = gcd(a, b) berechnen via erweitertem euklidischen Algorithmus
        2. Prüfen ob g | c
        3. Partikuläre Lösung: x0 = (c/g)*x', y0 = (c/g)*y'
           wobei a*x' + b*y' = g (aus Bezout-Identität)
        4. Allgemeine Lösung: x = x0 + (b/g)*t, y = y0 - (a/g)*t  für alle t ∈ ℤ

        Beispiel: 3x + 5y = 1
        - gcd(3,5) = 1, 1 | 1 ✓
        - Bezout: 3*2 + 5*(-1) = 1 → x0=2, y0=-1
        - Allgemein: x=2+5t, y=-1-3t

    @param a Koeffizient von x.
    @param b Koeffizient von y.
    @param c Rechte Seite.
    @return Tupel (x0, y0, gcd_ab) wobei (x0 + b/g*t, y0 - a/g*t) alle Lösungen sind,
            oder None wenn keine Lösung existiert.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # Erweiterter Euklidischer Algorithmus liefert g, x', y' mit a*x' + b*y' = g
    g, x_bezout, y_bezout = extended_gcd(a, b)

    # Notwendige und hinreichende Bedingung: g | c
    if c % g != 0:
        return None

    # Skalierungsfaktor für partikuläre Lösung
    factor = c // g

    # Partikuläre Lösung berechnen
    x0 = factor * x_bezout
    y0 = factor * y_bezout

    return (x0, y0, g)


def solve_quadratic_diophantine_pell(D: int, n_solutions: int = 5) -> list[tuple[int, int]]:
    """
    @brief Löst die Pell-Gleichung x² - D·y² = 1 (D kein perfektes Quadrat).
    @description
        Die Pell-Gleichung x² - D·y² = 1 hat unendlich viele ganzzahlige Lösungen,
        falls D kein perfektes Quadrat ist.

        Fundamentallösung (x₁, y₁) via Kettenbruchentwicklung von √D:
        - √D = a₀ + 1/(a₁ + 1/(a₂ + ...)) ist eine periodische Kettenbruchentwicklung
        - Die Fundamentallösung ergibt sich aus dem (r-1)-ten Konvergenten,
          wobei r die Periode des Kettenbruchs ist
        - Alle weiteren Lösungen via Rekursion:
          x_{k+1} = x₁·x_k + D·y₁·y_k
          y_{k+1} = x₁·y_k + y₁·x_k

        Historisch: Brahmagupta (7. Jh.), Bhaskara II (12. Jh.), Lagrange (1768)

    @param D Ganzzahl, kein perfektes Quadrat.
    @param n_solutions Anzahl der zu berechnenden Lösungen.
    @return Liste von (x, y) Tupeln, erste Lösungen der Pell-Gleichung.
    @raises ValueError Wenn D ein perfektes Quadrat ist.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # D darf kein perfektes Quadrat sein
    sqrt_D = int(math.isqrt(D))
    if sqrt_D * sqrt_D == D:
        raise ValueError(f"D={D} ist ein perfektes Quadrat, keine Pell-Gleichung")

    # Kettenbruchentwicklung von √D finden
    # a0 = floor(√D), dann periodisch
    a0 = sqrt_D

    # Konvergenten berechnen bis Pell-Gleichung erfüllt ist
    # Algorithmus: m, d, a initialisieren, dann iterieren
    m = 0
    d = 1
    a = a0

    # Konvergenten-Numerator und -Denominator berechnen
    p_prev, p_curr = 1, a0  # h_{-1}=1, h_0=a0
    q_prev, q_curr = 0, 1   # k_{-1}=0, k_0=1

    # Fundamentallösung durch Iteration finden
    for _ in range(10000):
        # Nächstes Kettenbruch-Glied berechnen
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d

        # Konvergenten aktualisieren
        p_prev, p_curr = p_curr, a * p_curr + p_prev
        q_prev, q_curr = q_curr, a * q_curr + q_prev

        # Pell-Gleichung prüfen: x²-D·y²=1?
        if p_curr * p_curr - D * q_curr * q_curr == 1:
            x1, y1 = p_curr, q_curr
            break
    else:
        raise RuntimeError(f"Fundamentallösung für D={D} nicht gefunden")

    # Alle Lösungen via Rekursionsformel erzeugen
    solutions = [(x1, y1)]
    xk, yk = x1, y1

    for _ in range(n_solutions - 1):
        # Rekursion: x_{k+1} = x1·x_k + D·y1·y_k, y_{k+1} = x1·y_k + y1·x_k
        xk_new = x1 * xk + D * y1 * yk
        yk_new = x1 * yk + y1 * xk
        xk, yk = xk_new, yk_new
        solutions.append((xk, yk))

    return solutions


def solve_pythagorean_triples(n: int) -> list[tuple[int, int, int]]:
    """
    @brief Findet alle primitiven Pythagoräischen Tripel (a, b, c) mit c ≤ n.
    @description
        Ein Pythagoräisches Tripel (a, b, c) erfüllt a² + b² = c².
        Es ist primitiv wenn gcd(a, b, c) = 1.

        Parametrische Darstellung (Euklid's Formel):
        Für m > k > 0, gcd(m,k) = 1 und m - k ungerade:
            a = m² - k²
            b = 2mk
            c = m² + k²

        Dies erzeugt alle primitiven Tripel genau einmal.

        Beispiele: (3,4,5), (5,12,13), (8,15,17), (7,24,25)

    @param n Obere Schranke für die Hypotenuse c.
    @return Sortierte Liste von (a, b, c) Tupeln mit a < b < c.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Beispiele:
    >>> solve_pythagorean_triples(20)
    [(3, 4, 5), (5, 12, 13), (8, 15, 17)]
    >>> solve_pythagorean_triples(5)
    [(3, 4, 5)]
    """
    triples = []

    # m von 2 bis ca. sqrt(n) iterieren (c = m²+k² ≤ n → m ≤ sqrt(n))
    m_max = int(math.isqrt(n)) + 1

    for m in range(2, m_max + 1):
        for k in range(1, m):
            # Bedingungen für primitive Tripel:
            # 1. gcd(m, k) = 1 (teilerfremd)
            # 2. m - k ist ungerade (verschiedene Parität)
            if gcd(m, k) != 1:
                continue
            if (m - k) % 2 == 0:
                continue

            # Euklid's Formel anwenden
            a = m * m - k * k
            b = 2 * m * k
            c = m * m + k * k

            # Schranke prüfen
            if c > n:
                break

            # a < b sicherstellen (sortiert)
            if a > b:
                a, b = b, a

            triples.append((a, b, c))

    # Nach c, dann a sortieren
    triples.sort(key=lambda t: (t[2], t[0]))
    return triples


def solve_diophantine_two_squares(n: int) -> list[tuple[int, int]]:
    """
    @brief Schreibt n als Summe zweier Quadrate: n = a² + b².
    @description
        Fermats Theorem der zwei Quadrate:
        Eine natürliche Zahl n ist genau dann als Summe zweier Quadrate darstellbar,
        wenn in der Primfaktorzerlegung jede Primzahl der Form 4k+3 mit gerader
        Vielfachheit vorkommt.

        Primzahlen p = 4k+1 sind immer als Summe zweier Quadrate darstellbar (Fermat).
        Primzahlen p = 4k+3 können nur in gerader Potenz auftreten.

        Algorithmus: Direktes Durchsuchen aller a von 0 bis floor(√n).
        Für jeden a prüfen ob n - a² ein perfektes Quadrat ist.

        Beispiele:
        - 25 = 0² + 5² = 3² + 4²
        - 50 = 1² + 7² = 5² + 5²
        - 3 → nicht darstellbar

    @param n Die zu zerlegende natürliche Zahl.
    @return Liste von (a, b) Tupeln mit a ≤ b und a² + b² = n, leer wenn nicht möglich.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    representations = []

    # a von 0 bis floor(√n) durchsuchen
    a_max = int(math.isqrt(n))

    for a in range(a_max + 1):
        remainder = n - a * a
        if remainder < 0:
            break

        # Prüfen ob remainder ein perfektes Quadrat ist
        b = int(math.isqrt(remainder))
        if b * b == remainder and a <= b:
            representations.append((a, b))

    return representations


def markov_numbers(n_terms: int = 20) -> list[int]:
    """
    @brief Berechnet Markov-Zahlen (Lösungen der Markov-Gleichung x² + y² + z² = 3xyz).
    @description
        Die Markov-Gleichung x² + y² + z² = 3xyz hat unendlich viele ganzzahlige
        Lösungen, die Markov-Tripel genannt werden. Die zugehörigen Zahlen heißen
        Markov-Zahlen.

        Baumstruktur der Lösungen:
        - Fundamentallösung: (1, 1, 1)
        - Aus jedem Tripel (x, y, z) erhält man durch "Spiegelung" ein neues:
          z' = 3xy - z  (wenn x² + y² + z² = 3xyz, dann auch mit z' statt z)
        - Aus (1,1,1): → (1,1,2) → (1,2,5) → (2,5,29) → ...

        Markov-Zahlen: 1, 2, 5, 13, 29, 34, 89, 169, 194, 233, ...

        Die Markov-Eindeutigkeitsvermutung (offen seit 1879) besagt:
        Jede Markov-Zahl bestimmt ihr Tripel eindeutig.

    @param n_terms Anzahl der zu berechnenden Markov-Zahlen.
    @return Sortierte Liste der ersten n_terms Markov-Zahlen.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    # BFS-Ansatz: Alle Tripel per Breitensuche erzeugen
    from collections import deque

    markov_set = set()
    # Starttripel
    queue = deque([(1, 1, 1)])
    visited = set()
    visited.add((1, 1, 1))

    while queue and len(markov_set) < n_terms * 3:
        x, y, z = queue.popleft()

        # Alle drei Zahlen des Tripels zur Menge hinzufügen
        markov_set.add(x)
        markov_set.add(y)
        markov_set.add(z)

        # Neue Tripel durch Spiegelung erzeugen
        # Spiegeln von z: z' = 3xy - z
        candidates = [
            tuple(sorted((3 * y * z - x, y, z))),
            tuple(sorted((x, 3 * x * z - y, z))),
            tuple(sorted((x, y, 3 * x * y - z))),
        ]

        for candidate in candidates:
            if candidate not in visited:
                visited.add(candidate)
                queue.append(candidate)

    # Sortiert zurückgeben
    result = sorted(markov_set)
    return result[:n_terms]


# ===========================================================================
# QUADRATISCHES REZIPROZITÄTSGESETZ
# ===========================================================================

def is_quadratic_residue(a: int, p: int) -> bool:
    """
    @brief Prüft ob a ein quadratischer Rest modulo der Primzahl p ist.
    @description
        a ist ein quadratischer Rest mod p wenn es ein x gibt mit x² ≡ a (mod p).

        Euler-Kriterium: a ist QR mod p ⟺ a^((p-1)/2) ≡ 1 (mod p)
        (für p ungerade Primzahl und p ∤ a)

        Das Euler-Kriterium folgt aus dem Kleinen Fermatschen Satz:
        a^(p-1) ≡ 1 (mod p) → (a^((p-1)/2))² ≡ 1 → a^((p-1)/2) ≡ ±1

        Falls a ≡ 0 (mod p): kein QR (Trivialfall)

    @param a Die zu prüfende Zahl.
    @param p Ungerade Primzahl.
    @return True wenn a ein quadratischer Rest mod p ist.
    @author Kurt Ingwer
    @lastModified 2026-03-10

    Beispiele:
    >>> is_quadratic_residue(2, 7)
    True
    >>> is_quadratic_residue(3, 7)
    False
    >>> is_quadratic_residue(1, 5)
    True
    """
    a = a % p
    # 0 ist kein quadratischer Rest im üblichen Sinne
    if a == 0:
        return False
    # Euler-Kriterium: a^((p-1)/2) ≡ 1 (mod p)
    return pow(a, (p - 1) // 2, p) == 1


def quadratic_residues(p: int) -> list[int]:
    """
    @brief Gibt alle quadratischen Reste modulo p zurück.
    @description
        Die quadratischen Reste mod p sind die Elemente von {1,...,p-1}
        die ein Quadrat modulo p sind.

        Eigenschaft: Es gibt genau (p-1)/2 quadratische Reste mod p (für Primzahl p).
        Diese sind die Werte {1², 2², ..., ((p-1)/2)²} mod p (alle verschieden).

        Beispiel p=7: 1²=1, 2²=4, 3²=2 mod 7 → QR={1,2,4}

    @param p Ungerade Primzahl.
    @return Sortierte Liste der quadratischen Reste mod p.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    residues = set()
    # Nur bis (p-1)/2 iterieren nötig, da k² ≡ (p-k)² mod p
    for i in range(1, p):
        residues.add((i * i) % p)
    # 0 entfernen falls vorhanden
    residues.discard(0)
    return sorted(residues)


def quadratic_reciprocity(p: int, q: int) -> dict[str, int | bool]:
    """
    @brief Wendet das Quadratische Reziprozitätsgesetz an.
    @description
        Das Quadratische Reziprozitätsgesetz (Gauß, 1796) beschreibt die Beziehung
        zwischen dem Legendre-Symbol (p|q) und (q|p):

            (p/q) · (q/p) = (-1)^{((p-1)/2) · ((q-1)/2)}

        Das bedeutet:
        - Wenn p ≡ 1 (mod 4) oder q ≡ 1 (mod 4): (p/q) = (q/p)
        - Wenn p ≡ q ≡ 3 (mod 4): (p/q) = -(q/p)

        Das Legendre-Symbol (a/p) = 1 wenn a QR mod p, -1 wenn NQR mod p, 0 wenn p|a.

    @param p Ungerade Primzahl.
    @param q Ungerade Primzahl (p ≠ q).
    @return Dict mit Legendre-Symbolen, Produkt, erwartetem Vorzeichen und Verifikation.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    def legendre_symbol(a: int, prime: int) -> int:
        """Berechnet das Legendre-Symbol (a/prime) via Euler-Kriterium."""
        a = a % prime
        if a == 0:
            return 0
        val = pow(a, (prime - 1) // 2, prime)
        # val ist 1 oder prime-1 (≡ -1)
        return 1 if val == 1 else -1

    # Legendre-Symbole berechnen
    lpq = legendre_symbol(p, q)  # (p/q)
    lqp = legendre_symbol(q, p)  # (q/p)

    # Vorzeichen gemäß Reziprozitätsgesetz
    exponent = ((p - 1) // 2) * ((q - 1) // 2)
    expected_sign = (-1) ** exponent

    # Produkt der beiden Legendre-Symbole
    product = lpq * lqp

    # Gesetz verifizieren
    law_verified = (product == expected_sign)

    return {
        'legendre_p_q': lpq,
        'legendre_q_p': lqp,
        'product': product,
        'expected_sign': expected_sign,
        'law_verified': law_verified
    }


def tonelli_shanks(n: int, p: int) -> int | None:
    """
    @brief Berechnet √n (mod p) via Tonelli-Shanks-Algorithmus.
    @description
        Der Tonelli-Shanks-Algorithmus findet x mit x² ≡ n (mod p).
        Er funktioniert für beliebige ungerade Primzahlen p.

        Spezialfall: p ≡ 3 (mod 4): x = n^((p+1)/4) mod p (direkte Formel)

        Allgemeiner Fall (Tonelli-Shanks):
        1. Schreibe p-1 = Q · 2^S mit Q ungerade
        2. Finde z mit z QNR mod p (quadratischer Nicht-Rest)
        3. Initialisiere M=S, c=z^Q, t=n^Q, R=n^((Q+1)/2) mod p
        4. Iteriere: Reduziere M bis Konvergenz

    @param n Zahl deren Quadratwurzel gesucht wird.
    @param p Ungerade Primzahl.
    @return x mit x² ≡ n (mod p), oder None wenn n kein QR mod p ist.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    n = n % p

    # Triviale Fälle
    if n == 0:
        return 0
    if p == 2:
        return n

    # Kein quadratischer Rest?
    if pow(n, (p - 1) // 2, p) != 1:
        return None

    # Spezialfall: p ≡ 3 (mod 4) → direkte Formel
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # Allgemeiner Tonelli-Shanks
    # Schritt 1: p-1 = Q * 2^S zerlegen
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    # Schritt 2: Quadratischen Nicht-Rest z finden
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    # Schritt 3: Initialisierung
    M = S
    c = pow(z, Q, p)       # c = z^Q mod p
    t = pow(n, Q, p)       # t = n^Q mod p
    R = pow(n, (Q + 1) // 2, p)  # R = n^((Q+1)/2) mod p

    # Schritt 4: Hauptschleife
    while True:
        if t == 0:
            return 0
        if t == 1:
            return R

        # i finden: kleinste i > 0 mit t^(2^i) ≡ 1 (mod p)
        i = 1
        temp = (t * t) % p
        while temp != 1:
            temp = (temp * temp) % p
            i += 1

        # Update
        b = pow(c, pow(2, M - i - 1), p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p


def cipolla_algorithm(n: int, p: int) -> int | None:
    """
    @brief Berechnet √n (mod p) via Cipolla-Algorithmus.
    @description
        Der Cipolla-Algorithmus (Michele Cipolla, 1903) ist eine Alternative
        zu Tonelli-Shanks zur Berechnung von Quadratwurzeln modulo Primzahlen.

        Algorithmus:
        1. Finde a, so dass a² - n ein QNR (quadratischer Nicht-Rest) mod p ist
        2. Arbeite in GF(p²) = GF(p)[x]/(x² - (a²-n))
           (erweitertes Körper über GF(p) mit irreduziblem Element ω mit ω² = a²-n)
        3. Berechne (a + ω)^((p+1)/2) in GF(p²)
        4. Das Ergebnis hat imaginären Teil 0 und reellen Teil = √n mod p

    @param n Zahl deren Quadratwurzel gesucht wird.
    @param p Ungerade Primzahl.
    @return x mit x² ≡ n (mod p), oder None wenn n kein QR mod p ist.
    @author Kurt Ingwer
    @lastModified 2026-03-10
    """
    n = n % p

    # Trivialfälle
    if n == 0:
        return 0
    if p == 2:
        return n

    # n muss quadratischer Rest sein
    if pow(n, (p - 1) // 2, p) != 1:
        return None

    # Schritt 1: a finden so dass a²-n ein QNR ist
    # (a²-n ist QNR ⟺ (a²-n)^((p-1)/2) ≡ -1 ≡ p-1 mod p)
    a = 0
    omega_sq = 0  # ω² = a² - n
    for a in range(p):
        omega_sq = (a * a - n) % p
        if omega_sq == 0:
            # a² ≡ n mod p → a ist direkt die Wurzel
            return a
        if pow(omega_sq, (p - 1) // 2, p) == p - 1:
            # omega_sq ist QNR → gefunden
            break
    else:
        return None  # Sollte nicht vorkommen

    # Schritt 2+3: Potenzierung in GF(p²)
    # Elemente der Form (x, y) entsprechen x + y*ω in GF(p²)
    # Multiplikation: (x1+y1·ω)·(x2+y2·ω) = (x1x2 + y1y2·ω²) + (x1y2+x2y1)·ω

    def mul_gf(x1: int, y1: int, x2: int, y2: int) -> tuple[int, int]:
        """Multiplikation in GF(p²)."""
        real = (x1 * x2 + y1 * y2 * omega_sq) % p
        imag = (x1 * y2 + x2 * y1) % p
        return real, imag

    # Schnelle Potenzierung: (a + ω)^((p+1)/2) mod p
    exp = (p + 1) // 2
    result_real, result_imag = 1, 0  # neutrale Element = 1
    base_real, base_imag = a % p, 1  # Basis = a + 1*ω

    while exp > 0:
        if exp % 2 == 1:
            result_real, result_imag = mul_gf(result_real, result_imag, base_real, base_imag)
        base_real, base_imag = mul_gf(base_real, base_imag, base_real, base_imag)
        exp //= 2

    # Imaginärteil muss 0 sein (Ergebnis liegt in GF(p))
    # result_real ist die gesuchte Quadratwurzel
    return result_real
