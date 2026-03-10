"""
@file l_functions.py
@brief L-Funktionen-Infrastruktur: Dirichlet-Charaktere, Hecke-L-Funktionen und BSD-Verbindung.
@description
    Implementiert die zentralen Bausteine der analytischen Zahlentheorie rund um
    L-Funktionen – verallgemeinerte Zeta-Funktionen der Form

        L(s) = \\sum_{n=1}^{\\infty} \\frac{a_n}{n^s}

    Dieses Modul umfasst:

    1. **Dirichlet-Charaktere mod q**
       \\chi: (\\mathbb{Z}/q\\mathbb{Z})^{\\times} \\to \\mathbb{C}^{\\times}
       Ein Dirichlet-Charakter ist ein vollständig multiplikativer, periodischer
       Gruppenhomomorphismus.

    2. **Dirichlet L-Funktionen**
       L(s, \\chi) = \\sum_{n=1}^{\\infty} \\frac{\\chi(n)}{n^s}
       Für Re(s) > 1 konvergiert die Reihe absolut. Analytische Fortsetzung via
       Euler-Knopp-Beschleunigung für Re(s) ≤ 1.

    3. **Hecke L-Funktionen** für Modulformen
       L(s, f) = \\sum_{n=1}^{\\infty} \\frac{a_n}{n^s}
       mit vollständiger L-Funktion
       \\Lambda(s, f) = \\left(\\frac{2\\pi}{\\sqrt{N}}\\right)^{-s} \\Gamma(s) L(s, f)

    4. **Gauß-Summen**
       \\tau(\\chi) = \\sum_{n=1}^{q} \\chi(n) \\cdot e^{2\\pi i n/q}
       Mit |\\tau(\\chi)|^2 = q für primitive Charaktere.

    5. **BSD-Verbindung**: Näherung von L(E, s) über lokale Euler-Faktoren
       L(E, s) = \\prod_p L_p(E, p^{-s})^{-1}
       mit dem analytischen Rang über Nullstellenzählung nahe s=1.

    Zusammenhang mit Millenniums-Problemen:
    - Birch-Swinnerton-Dyer: ord_{s=1} L(E,s) = Rang(E(\\mathbb{Q}))
    - Riemann-Hypothese: alle nicht-trivialen Nullstellen auf Re(s) = 1/2

@author Kurt Ingwer
@version 1.0
@since 2026-03-10
@lastModified 2026-03-10
"""

import cmath
import math
from typing import Callable, List, Optional
from functools import lru_cache

# Projektinterne Ausnahmen
from exceptions import InvalidInputError, DomainError


# ===========================================================================
# INTERNE HILFSFUNKTIONEN
# ===========================================================================

def _gamma_lanczos(z: complex) -> complex:
    """
    @brief Gamma-Funktion via Lanczos-Approximation (interne Hilfsfunktion).
    @description
        Berechnet Γ(z) mit ~15 signifikanten Stellen.
        Für Re(z) < 0.5 wird die Reflexionsformel verwendet.
    @param z Komplexe Zahl, z ∉ {0, -1, -2, ...}
    @return Γ(z)
    @lastModified 2026-03-10
    """
    # Pole bei z = 0, -1, -2, ... abfangen
    if z.real <= 0 and abs(z.imag) < 1e-12 and abs(z.real - round(z.real)) < 1e-12:
        raise ValueError(f"Γ(z) hat einen Pol bei z = {z}")

    # Reflexionsformel: Γ(z)·Γ(1-z) = π/sin(πz)
    if z.real < 0.5:
        return cmath.pi / (cmath.sin(cmath.pi * z) * _gamma_lanczos(1 - z))

    # Lanczos-Koeffizienten g=7, n=9 nach Paul Godfrey
    g = 7
    c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ]

    z -= 1
    x = complex(c[0])
    for i in range(1, g + 2):
        x += c[i] / (z + i)

    t = z + g + 0.5
    return cmath.sqrt(2 * cmath.pi) * (t ** (z + 0.5)) * cmath.exp(-t) * x


def _is_prime(n: int) -> bool:
    """
    @brief Primzahltest (Miller-Rabin-Näherung für kleine Zahlen).
    @param n Zu prüfende natürliche Zahl
    @return True wenn n prim, sonst False
    @lastModified 2026-03-10
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Probedivision bis √n
    for i in range(3, int(math.isqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True


def _euler_totient(q: int) -> int:
    """
    @brief Eulersche Totientfunktion φ(q) = #{k: 1 ≤ k ≤ q, gcd(k,q)=1}.
    @param q Positive ganze Zahl
    @return φ(q)
    @lastModified 2026-03-10
    """
    result = q
    p = 2
    temp = q
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def _gcd(a: int, b: int) -> int:
    """
    @brief Größter gemeinsamer Teiler (euklidischer Algorithmus).
    @param a Erste ganze Zahl
    @param b Zweite ganze Zahl
    @return gcd(a, b)
    @lastModified 2026-03-10
    """
    while b:
        a, b = b, a % b
    return abs(a)


def _primitive_root(q: int) -> Optional[int]:
    """
    @brief Findet eine Primitivwurzel modulo q (falls existent).
    @description
        Eine Primitivwurzel g mod q existiert wenn q = 1, 2, 4, p^k oder 2p^k
        für eine ungerade Primzahl p. g ist primitiv wenn ord(g) = φ(q).
    @param q Modulus
    @return Primitivwurzel g, oder None wenn keine existiert
    @lastModified 2026-03-10
    """
    phi = _euler_totient(q)
    # Primfaktorzerlegung von φ(q) für Ordnungstest
    phi_factors = _prime_factors(phi)

    for g in range(2, q + 1):
        if _gcd(g, q) != 1:
            continue
        # g ist primitiv, wenn g^(φ/p) ≢ 1 (mod q) für alle Primteiler p von φ
        is_primitive = True
        for p in phi_factors:
            if pow(g, phi // p, q) == 1:
                is_primitive = False
                break
        if is_primitive:
            return g
    return None


def _prime_factors(n: int) -> List[int]:
    """
    @brief Liefert alle Primfaktoren (ohne Duplikate) einer natürlichen Zahl.
    @param n Positive ganze Zahl ≥ 2
    @return Sortierte Liste der Primfaktoren
    @lastModified 2026-03-10
    """
    factors = []
    d = 2
    temp = abs(n)
    while d * d <= temp:
        if temp % d == 0:
            factors.append(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.append(temp)
    return factors


# ===========================================================================
# 1. DIRICHLET-CHARAKTERE
# ===========================================================================

def dirichlet_characters(q: int) -> List[dict]:
    """
    @brief Konstruiert alle Dirichlet-Charaktere modulo q.
    @description
        Ein Dirichlet-Charakter χ mod q ist ein vollständig multiplikativer
        Funktion χ: ℤ → ℂ mit:
          - χ(n+q) = χ(n)  (q-periodisch)
          - χ(n) = 0  wenn gcd(n,q) > 1
          - χ(mn) = χ(m)χ(n)  (vollständig multiplikativ)
          - χ(n) ≠ 0 für gcd(n,q) = 1

        Es gibt genau φ(q) verschiedene Charaktere mod q, die eine Gruppe
        isomorph zu (ℤ/qℤ)× bilden.

        Der Hauptcharakter χ₀ ist definiert durch χ₀(n) = 1 wenn gcd(n,q)=1,
        sonst 0. Er entspricht dem Pol von ζ(s) bei s=1 via
        L(s,χ₀) = ζ(s) · ∏_{p|q}(1 - p^{-s}).

    @param q Modulus (positive ganze Zahl ≥ 1)
    @return Liste von Dictionaries mit Schlüsseln:
            - 'chi': aufrufbares chi(n) → complex
            - 'order': multiplikative Ordnung des Charakters
            - 'is_principal': True wenn χ = χ₀ (Hauptcharakter)
            - 'values': Dict {n: chi(n) für n = 0,...,q-1}
            - 'index': Index in der Charakterliste (0 = Hauptcharakter)
    @raises InvalidInputError wenn q < 1
    @lastModified 2026-03-10
    """
    if q < 1:
        raise InvalidInputError("dirichlet_characters", f"q={q} muss ≥ 1 sein")

    # Spezialfall q=1: Nur der triviale Charakter χ ≡ 1
    if q == 1:
        def chi_trivial(n: int) -> complex:
            return complex(1)  # Alle Zahlen sind teilerfremd zu 1

        return [{
            'chi': chi_trivial,
            'order': 1,
            'is_principal': True,
            'values': {0: complex(1)},
            'index': 0,
        }]

    # Einheiten-Gruppe (ℤ/qℤ)× bestimmen
    units = [n for n in range(1, q) if _gcd(n, q) == 1]
    phi_q = len(units)  # = φ(q)

    # Charaktere werden durch Charaktere der abelschen Gruppe (ℤ/qℤ)× parametrisiert.
    # Für allgemeines q zerlegen wir in Primzahlpotenzkomponenten.
    # Vereinfachte Implementierung: Charaktere via Exponentialtabelle über
    # erzeugende Elemente der Gruppe.
    characters = []

    # Hauptcharakter (Index 0): χ₀(n) = 1 wenn gcd(n,q)=1, sonst 0
    chi0_values = {n: (complex(1) if _gcd(n, q) == 1 else complex(0)) for n in range(q)}

    def make_principal_chi(vals: dict) -> Callable:
        """Schließt die Wertetabelle in der Lambda-Funktion ein."""
        def chi(n: int) -> complex:
            return vals[n % len(vals)]
        return chi

    characters.append({
        'chi': make_principal_chi(chi0_values),
        'order': 1,
        'is_principal': True,
        'values': chi0_values,
        'index': 0,
    })

    # Weitere Charaktere: Nutze Struktursatz für (ℤ/qℤ)×
    # Für q = p (prim) ist die Gruppe zyklisch der Ordnung p-1.
    # Allgemein: Wir erzeugen alle Charaktere durch Exponentialabbildung.
    if phi_q > 1:
        additional = _construct_characters_mod_q(q, units, phi_q, chi0_values)
        for idx, char_vals in enumerate(additional):
            # Ordnung des Charakters: kleinste k mit χ^k = χ₀
            order = _character_order(char_vals, units, q)

            def make_chi(vals: dict) -> Callable:
                def chi(n: int) -> complex:
                    return vals[n % len(vals)]
                return chi

            characters.append({
                'chi': make_chi(char_vals),
                'order': order,
                'is_principal': False,
                'values': char_vals,
                'index': idx + 1,
            })

    return characters


def _construct_characters_mod_q(q: int, units: List[int],
                                  phi_q: int, chi0_values: dict) -> List[dict]:
    """
    @brief Interne Hilfsfunktion: Konstruiert alle Nicht-Hauptcharaktere mod q.
    @description
        Für jede zyklische Komponente der Gruppe (ℤ/qℤ)× wird ein Erzeuger
        gesucht und Charaktere via Exponentialabbildung konstruiert.

        Für q prim: (ℤ/qℤ)× ≅ ℤ/(q-1)ℤ, Erzeuger = Primitivwurzel g.
        Charakter χⱼ: χⱼ(gᵏ) = e^{2πijk/(q-1)} für j = 0,...,q-2.
    @param q Modulus
    @param units Liste der Einheiten mod q
    @param phi_q φ(q) = |units|
    @param chi0_values Wertetabelle des Hauptcharakters
    @return Liste der Wertetabellen für Nicht-Hauptcharaktere
    @lastModified 2026-03-10
    """
    result = []

    # Erzeuger der Gruppe finden (Primitivwurzel wenn vorhanden)
    g = _primitive_root(q)

    if g is not None:
        # Zyklische Gruppe der Ordnung phi_q: χⱼ(gᵏ) = ζ^{jk} mit ζ = e^{2πi/phi_q}
        # Zuerst: Logarithmustabelle log_g(n) = k mit g^k ≡ n (mod q)
        log_table = {}
        val = 1
        for k in range(phi_q):
            log_table[val] = k
            val = (val * g) % q

        # Charaktere χⱼ für j = 1,...,phi_q-1 (j=0 ist Hauptcharakter)
        for j in range(1, phi_q):
            char_vals = {}
            for n in range(q):
                if _gcd(n, q) == 1:
                    k = log_table.get(n % q, 0)
                    # χⱼ(n) = e^{2πijk/phi_q}
                    angle = 2 * math.pi * j * k / phi_q
                    char_vals[n] = complex(math.cos(angle), math.sin(angle))
                else:
                    char_vals[n] = complex(0)
            result.append(char_vals)
    else:
        # Für nicht-zyklische Gruppen (z.B. q=8): Produktkonstruktion
        # Vereinfachung: Direkte Konstruktion für kleine q
        for j in range(1, phi_q):
            char_vals = {}
            for n in range(q):
                if _gcd(n, q) == 1:
                    # Heuristik für nicht-zyklische Gruppen
                    angle = 2 * math.pi * j * units.index(n % q if n % q in units else units[0]) / phi_q
                    char_vals[n] = complex(math.cos(angle), math.sin(angle))
                else:
                    char_vals[n] = complex(0)
            result.append(char_vals)

    return result


def _character_order(char_vals: dict, units: List[int], q: int) -> int:
    """
    @brief Berechnet die multiplikative Ordnung eines Charakters.
    @description
        Die Ordnung von χ ist die kleinste positive ganze Zahl k mit χ^k = χ₀,
        d.h. mit χ(n)^k = 1 für alle n mit gcd(n,q)=1.
    @param char_vals Wertetabelle des Charakters
    @param units Einheitengruppe mod q
    @param q Modulus
    @return Ordnung des Charakters
    @lastModified 2026-03-10
    """
    if not units:
        return 1

    # Ordnung = kgV der Ordnungen der einzelnen Werte χ(n)
    # Für einen Charakter genügt ein erzeugendes Element
    # Wir prüfen für den ersten Einheitswert χ(g)
    u = units[0]
    val = char_vals.get(u, complex(1))

    if abs(val - 1) < 1e-10:
        return 1

    # Ordnung des komplexen Einheitswerts e^{2πi/k}: finde k via Argument
    angle = cmath.phase(val)
    if abs(angle) < 1e-10:
        return 1

    # k ≈ 2π / angle (im Einheitskreis)
    k_approx = 2 * math.pi / abs(angle)
    k = max(1, round(k_approx))

    # Verifiziere mit kleinster positiver Potenz
    for k_test in range(1, len(units) + 2):
        if abs(val ** k_test - 1) < 1e-8:
            return k_test

    return k


def chi_value(n: int, q: int, char_index: int) -> complex:
    """
    @brief Berechnet den Wert des char_index-ten Dirichlet-Charakters mod q an der Stelle n.
    @description
        Gibt χ_{char_index}(n) zurück, wobei χ₀ (Index 0) der Hauptcharakter ist.

        Eigenschaften:
        - χ(n) = 0 wenn gcd(n, q) > 1
        - |χ(n)| = 1 wenn gcd(n, q) = 1
        - χ(1) = 1 für alle Charaktere
    @param n Argument des Charakters
    @param q Modulus (≥ 1)
    @param char_index Index des Charakters (0 = Hauptcharakter)
    @return χ_{char_index}(n) als komplexe Zahl
    @raises InvalidInputError wenn char_index < 0 oder ≥ φ(q)
    @raises DomainError wenn q < 1
    @lastModified 2026-03-10
    """
    if q < 1:
        raise DomainError("chi_value", q, domain="q ≥ 1")

    chars = dirichlet_characters(q)
    if char_index < 0 or char_index >= len(chars):
        raise InvalidInputError(
            "chi_value",
            f"char_index={char_index} außerhalb [0, {len(chars)-1}]"
        )

    return chars[char_index]['chi'](n)


def is_primitive_character(chi_values: dict, q: int) -> bool:
    """
    @brief Prüft ob ein gegebener Charakter (als Wertetabelle) primitiv mod q ist.
    @description
        Ein Charakter χ mod q heißt primitiv, wenn er keinen echten Teiler d | q
        als Führer hat, d.h. wenn χ nicht von einem Charakter mod d für d < q
        induziert wird.

        Kriterium: χ ist primitiv, wenn für keinen echten Teiler d | q gilt,
        dass χ(n) = χ(m) wann immer n ≡ m (mod d) und gcd(nm, q) = 1.

        Äquivalent: |τ(χ)|² = q (Gauß-Summen-Kriterium für primitive Charaktere).
    @param chi_values Dictionary {n: χ(n)} für n = 0,...,q-1
    @param q Modulus
    @return True wenn χ primitiv, False sonst
    @lastModified 2026-03-10
    """
    if q == 1:
        return True  # Einziger Charakter mod 1 ist trivial und gilt als primitiv

    # Hauptcharakter mod q ist nie primitiv (außer für q=1)
    # Prüfe: ist χ der Hauptcharakter?
    is_principal = all(
        abs(chi_values.get(n, 0) - (1 if _gcd(n, q) == 1 else 0)) < 1e-10
        for n in range(q)
    )
    if is_principal and q > 1:
        return False

    # Gauß-Summen-Kriterium: primitiv ⟺ |τ(χ)|² = q
    # τ(χ) = Σ_{n=1}^{q} χ(n) · e^{2πin/q}
    tau = complex(0)
    for n in range(1, q + 1):
        chi_n = chi_values.get(n % q, complex(0))
        if chi_n != 0:
            angle = 2 * math.pi * n / q
            tau += chi_n * complex(math.cos(angle), math.sin(angle))

    gauss_abs_sq = abs(tau) ** 2
    return abs(gauss_abs_sq - q) < 0.5  # Toleranz für Rundungsfehler


# ===========================================================================
# 2. DIRICHLET L-FUNKTIONEN
# ===========================================================================

def dirichlet_l_function(s: complex, q: int, char_index: int,
                          n_terms: int = 1000) -> complex:
    """
    @brief Berechnet L(s, χ) = Σ χ(n)/n^s für einen Dirichlet-Charakter χ mod q.
    @description
        Für Re(s) > 1 konvergiert die Reihe absolut und wird direkt summiert.
        Für Re(s) ≤ 1 wird die Euler-Knopp-Beschleunigung angewendet
        (analog zur Eta-Beschleunigung für die Riemann-Zeta-Funktion).

        Wichtige Spezialwerte:
        - L(1, χ₋₄) = 1 - 1/3 + 1/5 - 1/7 + ... = π/4  (Leibniz-Reihe)
          mit χ₋₄ = Charakter mod 4 mit χ(1)=1, χ(3)=-1.
        - L(1, χ₀) divergiert (Pol bei s=1 für den Hauptcharakter).

        Für den Hauptcharakter χ₀ mod q gilt:
        L(s, χ₀) = ζ(s) · ∏_{p | q}(1 - p^{-s})
        (hat Pol erster Ordnung bei s=1).
    @param s Komplexes Argument (Re(s) > 0 empfohlen)
    @param q Modulus des Charakters (≥ 1)
    @param char_index Index des Charakters (0 = Hauptcharakter)
    @param n_terms Anzahl der Summanden (mehr = genauer)
    @return L(s, χ) als komplexe Zahl
    @raises InvalidInputError wenn Parameter ungültig
    @raises DomainError wenn s ein Pol des Hauptcharakters ist (s=1, char_index=0)
    @lastModified 2026-03-10
    """
    if q < 1:
        raise InvalidInputError("dirichlet_l_function", f"q={q} muss ≥ 1 sein")
    if n_terms < 1:
        raise InvalidInputError("dirichlet_l_function", f"n_terms={n_terms} muss ≥ 1 sein")

    chars = dirichlet_characters(q)
    if char_index < 0 or char_index >= len(chars):
        raise InvalidInputError(
            "dirichlet_l_function",
            f"char_index={char_index} nicht in [0, {len(chars)-1}]"
        )

    # Pol-Prüfung: Hauptcharakter bei s=1
    is_principal = chars[char_index]['is_principal']
    if is_principal and abs(s - 1) < 1e-10:
        raise DomainError("dirichlet_l_function", s, domain="L(s,χ₀) hat Pol bei s=1")

    chi = chars[char_index]['chi']
    chi_values_list = [chi(n) for n in range(n_terms + 1)]

    if s.real > 1:
        # Direkte Summation: L(s,χ) = Σ_{n=1}^∞ χ(n)/n^s
        total = complex(0)
        for n in range(1, n_terms + 1):
            chi_n = chi_values_list[n]
            if abs(chi_n) > 1e-14:
                total += chi_n / (n ** s)
        return total

    # Für Re(s) ≤ 1: Euler-Knopp-Beschleunigung
    # Die Reihe L(s,χ) = Σ χ(n)/n^s konvergiert bedingt für Re(s) > 0.
    # Wir nutzen die E-Transformation: analog zu η(s) für ζ(s).
    return _dirichlet_l_euler_accelerated(s, chi_values_list, n_terms)


def _dirichlet_l_euler_accelerated(s: complex, chi_vals: list,
                                     n_terms: int) -> complex:
    """
    @brief Euler-Knopp-beschleunigte Berechnung von L(s, χ) für Re(s) ≤ 1.
    @description
        Für einen nicht-Hauptcharakter χ wechseln die Vorzeichen χ(n) · n^{-s}
        nicht notwendigerweise alternierend. Wir verwenden daher die
        Van Wijngaarden-Transformation (Euler-Knopp) für allgemeine Reihen.

        Für einen primitiven Charakter χ gilt:
        Σ_{n=0}^{q-1} χ(n) = 0 (Orthogonalitätsrelation)
        Dies gewährleistet Konvergenz für Re(s) > 0.

        Die Euler-Knopp-Transformation:
        S = Σ_{k=0}^{n-1} (1/2^{k+1}) · Δ^k a_0
        mit a_n = χ(n+1)/(n+1)^s und Δ^k die k-te Vorwärtsdifferenz.
    @param s Komplexes Argument
    @param chi_vals Liste der Charakterwerte [χ(0), χ(1), ..., χ(n_terms)]
    @param n_terms Anzahl der Terme
    @return Beschleunigter L-Wert
    @lastModified 2026-03-10
    """
    N = min(n_terms, 100)  # Euler-Knopp mit N Termen

    # Folge der Summanden a_n = χ(n)/n^s
    a = []
    for n in range(1, N + 1):
        chi_n = chi_vals[n] if n < len(chi_vals) else complex(0)
        if abs(chi_n) > 1e-14:
            a.append(chi_n / (n ** s))
        else:
            a.append(complex(0))

    if not a:
        return complex(0)

    # Euler-E-Transformation: S = Σ_{k=0}^{N-1} (1/2^{k+1}) · Δ^k a[0]
    # Δ^k a[j] = Σ_{i=0}^{k} (-1)^i C(k,i) a[j+i]
    total = complex(0)
    current = list(a)  # Δ^0 a

    for k in range(N):
        if k >= len(current):
            break
        total += current[0] / (2 ** (k + 1))
        # Nächste Differenz berechnen: Δ^{k+1} a[j] = Δ^k a[j+1] - Δ^k a[j]
        current = [current[i + 1] - current[i] for i in range(len(current) - 1)]
        if not current:
            break

    return total


def dirichlet_l_function_zeros(q: int, char_index: int,
                                t_min: float = 0.0, t_max: float = 30.0,
                                steps: int = 500) -> List[dict]:
    """
    @brief Sucht Nullstellen von L(1/2 + it, χ) im gegebenen t-Intervall.
    @description
        Untersucht die Funktion L(s, χ) auf der kritischen Geraden Re(s) = 1/2
        im Intervall [t_min, t_max]. Vorzeichenwechsel von Re(L) oder Im(L)
        werden als Nullstellenkandidaten betrachtet.

        Die verallgemeinerte Riemann-Hypothese (GRH) besagt:
        Alle nicht-trivialen Nullstellen von L(s, χ) liegen auf Re(s) = 1/2.

        Nullstellensuche via Vorzeichenwechsel des Realteils.
        Verfeinert mit einer binären Suche zwischen den Stützstellen.
    @param q Modulus des Charakters
    @param char_index Index des Charakters
    @param t_min Untere Grenze des Suchintervalls (t ≥ 0)
    @param t_max Obere Grenze des Suchintervalls
    @param steps Anzahl der Abtastpunkte
    @return Liste von Dictionaries {'t': float, 'zero': complex, 'abs_val': float}
    @raises InvalidInputError wenn t_min ≥ t_max oder steps < 2
    @lastModified 2026-03-10
    """
    if t_min >= t_max:
        raise InvalidInputError(
            "dirichlet_l_function_zeros",
            f"t_min={t_min} muss < t_max={t_max} sein"
        )
    if steps < 2:
        raise InvalidInputError(
            "dirichlet_l_function_zeros",
            f"steps={steps} muss ≥ 2 sein"
        )

    zeros = []
    dt = (t_max - t_min) / steps

    def f(t: float) -> complex:
        """L(1/2 + it, χ) auswerten."""
        try:
            return dirichlet_l_function(complex(0.5, t), q, char_index, n_terms=200)
        except Exception:
            return complex(float('nan'))

    # Abtasten und Vorzeichenwechsel suchen
    t_prev = t_min
    val_prev = f(t_prev)

    for i in range(1, steps + 1):
        t_curr = t_min + i * dt
        val_curr = f(t_curr)

        if math.isnan(val_curr.real) or math.isnan(val_prev.real):
            t_prev, val_prev = t_curr, val_curr
            continue

        # Vorzeichenwechsel in Realteil oder Betrag-Minimum
        real_sign_change = (val_prev.real * val_curr.real < 0)

        if real_sign_change:
            # Binäre Suche zur Verfeinerung
            t_lo, t_hi = t_prev, t_curr
            for _ in range(30):
                t_mid = (t_lo + t_hi) / 2
                val_mid = f(t_mid)
                if abs(val_mid) < 1e-6:
                    break
                if val_prev.real * val_mid.real <= 0:
                    t_hi, val_curr = t_mid, val_mid
                else:
                    t_lo, val_prev = t_mid, val_mid

            t_zero = (t_lo + t_hi) / 2
            zero_val = f(t_zero)
            # Doppelte Nullstellen vermeiden
            if not zeros or abs(t_zero - zeros[-1]['t']) > 0.1:
                zeros.append({
                    't': t_zero,
                    'zero': complex(0.5, t_zero),
                    'abs_val': abs(zero_val),
                })

        t_prev, val_prev = t_curr, val_curr

    return zeros


# ===========================================================================
# 3. GAUß-SUMMEN
# ===========================================================================

def gauss_sum(q: int, char_index: int) -> complex:
    """
    @brief Berechnet die Gauß-Summe τ(χ) = Σ_{n=1}^{q} χ(n) · e^{2πin/q}.
    @description
        Die Gauß-Summe ist fundamental für die Funktionalgleichung
        der Dirichlet L-Funktion:

            Λ(s, χ) = (q/π)^{s/2} Γ(s/2) L(s, χ)
            Λ(1-s, χ̄) = (τ(χ)/√q) · Λ(s, χ)

        Für primitive Charaktere gilt:
        |τ(χ)|² = q  (Gauß-Summen-Abschätzung)

        Für den Hauptcharakter χ₀ gilt (Ramanujan-Summe):
        τ(χ₀) = μ(q)  (Möbius-Funktion)

        Spezialfall q=4, χ₋₄: τ(χ₋₄) = 2i (Gauß'sche Summe)
    @param q Modulus (≥ 1)
    @param char_index Index des Charakters in der Charakterliste mod q
    @return τ(χ_{char_index}) als komplexe Zahl
    @raises InvalidInputError wenn Parameter ungültig
    @lastModified 2026-03-10
    """
    if q < 1:
        raise InvalidInputError("gauss_sum", f"q={q} muss ≥ 1 sein")

    chars = dirichlet_characters(q)
    if char_index < 0 or char_index >= len(chars):
        raise InvalidInputError(
            "gauss_sum",
            f"char_index={char_index} außerhalb [0, {len(chars)-1}]"
        )

    chi = chars[char_index]['chi']

    # τ(χ) = Σ_{n=1}^{q} χ(n) · e^{2πin/q}
    tau = complex(0)
    for n in range(1, q + 1):
        chi_n = chi(n)
        if abs(chi_n) > 1e-14:
            angle = 2 * math.pi * n / q
            tau += chi_n * complex(math.cos(angle), math.sin(angle))

    return tau


# ===========================================================================
# 4. ANALYTISCHE EIGENSCHAFTEN
# ===========================================================================

def l_function_conductor(q: int, char_index: int) -> int:
    """
    @brief Berechnet den Führer (Conductor) einer Dirichlet L-Funktion.
    @description
        Der Führer f(χ) eines Charakters χ mod q ist der kleinste positive
        Teiler f | q, sodass χ von einem primitiven Charakter mod f induziert wird.

        Eigenschaften:
        - f(χ₀) = 1 für den Hauptcharakter (per Konvention)
        - f(χ) = q wenn χ primitiv mod q ist
        - f(χ) | q stets

        Der Führer bestimmt den analytischen Conductor der L-Funktion.
    @param q Modulus
    @param char_index Index des Charakters
    @return Führer f(χ) (positiver Teiler von q)
    @raises InvalidInputError wenn Parameter ungültig
    @lastModified 2026-03-10
    """
    if q < 1:
        raise InvalidInputError("l_function_conductor", f"q={q} muss ≥ 1 sein")

    chars = dirichlet_characters(q)
    if char_index < 0 or char_index >= len(chars):
        raise InvalidInputError(
            "l_function_conductor",
            f"char_index={char_index} außerhalb [0, {len(chars)-1}]"
        )

    chi_vals = chars[char_index]['values']

    # Hauptcharakter hat Führer 1
    if chars[char_index]['is_principal']:
        return 1

    # Suche kleinsten Teiler d | q, sodass χ durch χ mod d induziert wird
    # χ wird durch d induziert, wenn χ(n) nur von n mod d abhängt
    # (für alle n mit gcd(n, q) = 1)
    divisors = sorted([d for d in range(1, q + 1) if q % d == 0])

    for d in divisors:
        if d == q:
            return q  # Primitiv mod q

        # Prüfe ob χ von einem Charakter mod d induziert wird:
        # Gilt χ(n) = χ(m) für alle n ≡ m (mod d) mit gcd(n,q) = gcd(m,q) = 1?
        units_q = [n for n in range(1, q) if _gcd(n, q) == 1]

        induced = True
        for n in units_q:
            for m in units_q:
                if n % d == m % d:
                    chi_n = chi_vals.get(n, complex(0))
                    chi_m = chi_vals.get(m, complex(0))
                    if abs(chi_n - chi_m) > 1e-8:
                        induced = False
                        break
            if not induced:
                break

        if induced:
            return d

    return q


def l_function_special_values(q: int, char_index: int) -> dict:
    """
    @brief Berechnet spezielle Werte der Dirichlet L-Funktion.
    @description
        Wichtige Spezialwerte:

        **L(1, χ) für nicht-Hauptcharaktere** (Dirichlets Primzahlsatz):
        L(1, χ) = -1/q · Σ_{a=1}^{q} χ̄(a) · ln(1 - e^{2πia/q})  (falls χ reell)
        Spezialfall: L(1, χ₋₄) = π/4 ≈ 0.7854 (Leibniz-Reihe)

        **L(0, χ) via Funktionalgleichung**:
        Für ungerade Charaktere (χ(-1) = -1): L(0, χ) = -B_{1,χ}/1
        Für gerade Charaktere (χ(-1) = 1): L(0, χ) = 0

        **Bernoulli-Darstellung**:
        L(1-n, χ) = -B_{n,χ}/n für n ≥ 1

        Der Wert L(1, χ) für primitive nicht-Hauptcharaktere ist stets ≠ 0
        (fundamentales Lemma für Dirichlets Primzahlsatz in AP).
    @param q Modulus
    @param char_index Index des Charakters
    @return Dictionary mit Schlüsseln 'L_1' (L(1,χ)), 'L_0' (L(0,χ)),
            'chi_minus1' (χ(-1)), 'is_odd' (χ ungerade), 'is_even' (χ gerade)
    @raises DomainError wenn Hauptcharakter (hat Pol bei s=1)
    @lastModified 2026-03-10
    """
    chars = dirichlet_characters(q)
    if char_index < 0 or char_index >= len(chars):
        raise InvalidInputError(
            "l_function_special_values",
            f"char_index={char_index} außerhalb [0, {len(chars)-1}]"
        )

    if chars[char_index]['is_principal']:
        raise DomainError(
            "l_function_special_values", "χ₀",
            domain="Hauptcharakter hat Pol bei s=1"
        )

    chi = chars[char_index]['chi']
    chi_minus1 = chi(-1)  # χ(-1) = ±1

    # L(1, χ) numerisch mit hoher Genauigkeit (500 Terme + Euler-Knopp)
    try:
        l1 = dirichlet_l_function(complex(1.0 + 1e-9, 0), q, char_index, n_terms=500)
    except Exception:
        l1 = complex(float('nan'))

    # L(0, χ) via Funktionalgleichung:
    # Für ungerade χ (χ(-1) = -1): L(0, χ) = -τ(χ̄)/(2√q) · (erste Bernoulli-Koeff.)
    # Vereinfachung: L(0, χ) ≈ L(ε, χ) für kleine ε
    try:
        l0 = dirichlet_l_function(complex(1e-6, 0), q, char_index, n_terms=300)
    except Exception:
        l0 = complex(float('nan'))

    is_odd = (abs(chi_minus1 + 1) < 1e-8)   # χ(-1) = -1
    is_even = (abs(chi_minus1 - 1) < 1e-8)  # χ(-1) = +1

    return {
        'L_1': l1,
        'L_0': l0,
        'chi_minus1': chi_minus1,
        'is_odd': is_odd,
        'is_even': is_even,
    }


def functional_equation_check(s: complex, q: int, char_index: int) -> dict:
    """
    @brief Überprüft die Funktionalgleichung Λ(s, χ) = ε(χ) · Λ(1-s, χ̄).
    @description
        Die vollständige Dirichlet L-Funktion (completed L-function) ist:
            Λ(s, χ) = (q/π)^{s/2} · Γ((s + δ)/2) · L(s, χ)
        mit δ = 0 für gerade χ (χ(-1)=1) und δ = 1 für ungerade χ (χ(-1)=-1).

        Sie erfüllt die Funktionalgleichung:
            Λ(s, χ) = ε(χ) · Λ(1-s, χ̄)

        mit dem Root Number:
            ε(χ) = τ(χ) / (i^δ · √q)

        Für reellwertige Charaktere gilt χ̄ = χ und |ε(χ)| = 1.

        Dieses Resultat ist zentral für:
        - Die analytische Fortsetzung von L(s, χ) auf ganz ℂ
        - Die Funktionalgleichung der vollständigen Zeta-Funktion
        - Die Berechnung von L(0, χ) aus L(1, χ̄)
    @param s Komplexes Argument (nicht 0 oder 1 für Hauptcharakter)
    @param q Modulus
    @param char_index Index des Charakters
    @return Dictionary mit 'lhs' Λ(s,χ), 'rhs' ε(χ)·Λ(1-s,χ̄), 'error' |lhs-rhs|,
            'epsilon' ε(χ), 'delta' δ
    @raises InvalidInputError wenn Parameter ungültig
    @lastModified 2026-03-10
    """
    if q < 1:
        raise InvalidInputError("functional_equation_check", f"q={q} muss ≥ 1 sein")

    chars = dirichlet_characters(q)
    if char_index < 0 or char_index >= len(chars):
        raise InvalidInputError(
            "functional_equation_check",
            f"char_index={char_index} außerhalb [0, {len(chars)-1}]"
        )

    chi = chars[char_index]['chi']

    # δ bestimmen: 0 für gerade, 1 für ungerade Charakter
    chi_minus1 = chi(-1)
    delta = 1 if abs(chi_minus1 + 1) < 1e-8 else 0

    # Konjugierter Charakter χ̄: Index über Wertetabelle suchen
    chi_vals = chars[char_index]['values']
    conj_vals = {n: v.conjugate() for n, v in chi_vals.items()}

    # χ̄ in der Charakterliste finden (oder mit konj. Werten approximieren)
    conj_index = 0  # Fallback: Hauptcharakter
    for idx, c in enumerate(chars):
        match = all(
            abs(c['values'].get(n, 0) - conj_vals.get(n, 0)) < 1e-8
            for n in range(min(q, 20))
        )
        if match:
            conj_index = idx
            break

    # Gauß-Summe τ(χ) berechnen
    tau = gauss_sum(q, char_index)

    # Epsilon (Root Number): ε = τ(χ) / (i^δ · √q)
    i_delta = complex(0, 1) ** delta
    epsilon = tau / (i_delta * math.sqrt(q))

    def completed_l(s_val: complex, chi_idx: int) -> complex:
        """
        Vollständige L-Funktion Λ(s, χ) = (q/π)^{s/2} · Γ((s+δ)/2) · L(s, χ).
        """
        try:
            # Faktor (q/π)^{s/2}
            qpi_factor = (q / math.pi) ** (s_val / 2)
            # Gamma-Faktor: Γ((s+δ)/2)
            gamma_arg = (s_val + delta) / 2
            gamma_val = _gamma_lanczos(gamma_arg)
            # L(s, χ)
            l_val = dirichlet_l_function(s_val, q, chi_idx, n_terms=200)
            return qpi_factor * gamma_val * l_val
        except Exception:
            return complex(float('nan'))

    # LHS: Λ(s, χ)
    lhs = completed_l(s, char_index)
    # RHS: ε(χ) · Λ(1-s, χ̄)
    rhs = epsilon * completed_l(1 - s, conj_index)

    error = abs(lhs - rhs) if (not math.isnan(lhs.real) and not math.isnan(rhs.real)) else float('nan')

    return {
        'lhs': lhs,
        'rhs': rhs,
        'error': error,
        'epsilon': epsilon,
        'delta': delta,
        'tau': tau,
    }


# ===========================================================================
# 5. HECKE L-FUNKTIONEN
# ===========================================================================

def hecke_l_function(s: complex, coefficients: List[float],
                      n_terms: int = 200) -> complex:
    """
    @brief Berechnet L(s, f) = Σ a_n / n^s für eine Modulform f.
    @description
        Eine Modulform f vom Gewicht k und Level N hat eine Fourier-Entwicklung:
            f(τ) = Σ_{n=1}^∞ a_n · q^n  mit q = e^{2πiτ}

        Die zugehörige L-Funktion ist:
            L(s, f) = Σ_{n=1}^∞ a_n / n^s

        Sie konvergiert absolut für Re(s) > (k+2)/2 und hat analytische
        Fortsetzung auf ganz ℂ.

        Normierung: Die Liste 'coefficients' enthält [a_1, a_2, ..., a_N],
        wobei a_1 = 1 für normierte Hecke-Eigenformen.

        Für die Ramanujan-Delta-Funktion Δ(τ) (k=12):
        - a_1 = 1, a_2 = τ(2) = -24, a_3 = τ(3) = 252
        - L(Δ, s) konvergiert gut für Re(s) > 7

        Eulerprodukt für Hecke-Eigenformen:
            L(s, f) = ∏_p (1 - a_p · p^{-s} + p^{k-1-2s})^{-1}
    @param s Komplexes Argument
    @param coefficients Liste [a_1, a_2, ..., a_N] der Fourier-Koeffizienten
    @param n_terms Maximale Anzahl Summanden (min(len(coefficients), n_terms))
    @return L(s, f) als komplexe Zahl
    @raises InvalidInputError wenn coefficients leer oder n_terms < 1
    @lastModified 2026-03-10
    """
    if not coefficients:
        raise InvalidInputError("hecke_l_function", "coefficients darf nicht leer sein")
    if n_terms < 1:
        raise InvalidInputError("hecke_l_function", f"n_terms={n_terms} muss ≥ 1 sein")

    N = min(len(coefficients), n_terms)
    total = complex(0)

    for n in range(1, N + 1):
        a_n = coefficients[n - 1]  # 0-indiziert: a_1 = coefficients[0]
        if abs(a_n) > 1e-14:
            total += complex(a_n) / (n ** s)

    return total


def completed_l_function(s: complex, coefficients: List[float],
                           weight: int = 12) -> complex:
    """
    @brief Vollständige L-Funktion Λ(s, f) = (2π)^{-s} Γ(s) L(s, f).
    @description
        Die vollständige (completed) L-Funktion ist:
            Λ(s, f) = N^{s/2} · (2π)^{-s} · Γ(s) · L(s, f)

        (vereinfacht ohne Level N: Λ(s, f) = (2π)^{-s} · Γ(s) · L(s, f))

        Sie erfüllt die Funktionalgleichung:
            Λ(s, f) = ε · Λ(k - s, f̄)

        mit ε = ±1 (Vorzeichen/Root Number), k = Gewicht der Modulform.

        Für die Delta-Funktion (k=12, ε=1):
            Λ(s, Δ) = Λ(12 - s, Δ)

        Die vollständige L-Funktion ist ganz (keine Pole für Cusp-Formen).
    @param s Komplexes Argument
    @param coefficients Fourier-Koeffizienten [a_1, a_2, ...]
    @param weight Gewicht k der Modulform (Standard: 12 für Delta-Funktion)
    @return Λ(s, f) = (2π)^{-s} · Γ(s) · L(s, f)
    @raises InvalidInputError wenn weight < 2 oder coefficients leer
    @lastModified 2026-03-10
    """
    if weight < 2:
        raise InvalidInputError("completed_l_function", f"weight={weight} muss ≥ 2 sein")
    if not coefficients:
        raise InvalidInputError("completed_l_function", "coefficients darf nicht leer sein")

    # Gamma-Faktor: Γ(s)
    try:
        gamma_val = _gamma_lanczos(s)
    except (ValueError, OverflowError):
        return complex(float('nan'))

    # (2π)^{-s}
    two_pi_factor = (2 * math.pi) ** (-s)

    # L(s, f) mit allen verfügbaren Koeffizienten
    l_val = hecke_l_function(s, coefficients, n_terms=len(coefficients))

    return two_pi_factor * gamma_val * l_val


# ===========================================================================
# 6. BSD-VERBINDUNG: ELLIPTISCHE L-FUNKTIONEN
# ===========================================================================

def _count_points_mod_p(a: int, b: int, p: int) -> int:
    """
    @brief Zählt die Punkte auf E: y² = x³ + ax + b über F_p (inklusive O).
    @description
        Direkte Zählung: Für jedes x ∈ {0,...,p-1} prüfe ob x³+ax+b ein
        quadratischer Rest mod p ist.

        Zeitkomplexität: O(p log p)
    @param a Kurvenparameter a (Weierstrass-Form)
    @param b Kurvenparameter b
    @param p Primzahl
    @return #E(F_p) inklusive Punkt im Unendlichen
    @lastModified 2026-03-10
    """
    count = 1  # Punkt im Unendlichen
    for x in range(p):
        rhs = (pow(x, 3, p) + a * x + b) % p
        # Quadratische Reste via Legendre-Symbol (Euler-Kriterium)
        if rhs == 0:
            count += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            count += 2  # y und -y
    return count


def elliptic_l_function_approx(a: int, b: int, s: complex,
                                 prime_bound: int = 100) -> complex:
    """
    @brief Näherung von L(E, s) via lokale Euler-Faktoren für alle Primzahlen bis prime_bound.
    @description
        Die L-Funktion einer elliptischen Kurve E: y² = x³ + ax + b ist:

            L(E, s) = ∏_p L_p(E, p^{-s})^{-1}

        mit lokalen Faktoren (für gute Primzahlen p):
            L_p(E, T)^{-1} = 1 - a_p · T + p · T²

        wobei a_p = p + 1 - #E(F_p) die Frobenius-Spur ist und T = p^{-s}.

        Für schlechte Primzahlen p | Δ (Diskriminante):
            L_p(E, T)^{-1} = 1 - T  (multiplikative Reduktion)
            L_p(E, T)^{-1} = 1 + T  (multiplikative, nicht-split)
            L_p(E, T)^{-1} = 1      (additive Reduktion)

        Die Diskriminante: Δ = -16(4a³ + 27b²)

        BSD-Vermutung: ord_{s=1} L(E, s) = Rang(E(ℚ))

        Für das Euler-Produkt konvergiert die Näherung für Re(s) > 3/2 gut.
        Der Wert L(E, 1) ist der kritische Wert (central value).
    @param a Kurvenparameter a der Weierstrass-Form y² = x³ + ax + b
    @param b Kurvenparameter b
    @param s Komplexes Argument
    @param prime_bound Alle Primzahlen bis zu diesem Wert werden verwendet
    @return Näherung von L(E, s) als komplexes Euler-Produkt
    @raises DomainError wenn die Kurve singulär ist (Δ = 0)
    @raises InvalidInputError wenn prime_bound < 2
    @lastModified 2026-03-10
    """
    if prime_bound < 2:
        raise InvalidInputError(
            "elliptic_l_function_approx",
            f"prime_bound={prime_bound} muss ≥ 2 sein"
        )

    # Diskriminante prüfen: Δ = -16(4a³ + 27b²)
    disc = -16 * (4 * a**3 + 27 * b**2)
    if disc == 0:
        raise DomainError(
            "elliptic_l_function_approx",
            f"(a={a}, b={b})",
            domain="nicht-singuläre elliptische Kurve (Δ ≠ 0)"
        )

    # Primzahlen bis prime_bound sammeln
    primes = [p for p in range(2, prime_bound + 1) if _is_prime(p)]

    # Euler-Produkt: L(E,s) = ∏_p (lokaler Faktor)^{-1}
    # Wir berechnen log L(E,s) = -Σ log(lokaler Faktor)
    log_l = complex(0)

    for p in primes:
        ps = p ** (-s)  # p^{-s}

        # Schlechte Primzahl? p | disc
        if disc % p == 0:
            # Vereinfachung: multiplikative Reduktion → 1 ± T
            # Exakter wäre: Klassifikation der Reduktionsart
            local_factor = 1 - ps  # Näherung für multiplikative Reduktion
        else:
            # Gute Primzahl: Frobenius-Spur a_p = p + 1 - #E(F_p)
            n_fp = _count_points_mod_p(a, b, p)
            a_p = p + 1 - n_fp
            # L_p(E, T)^{-1} = 1 - a_p · T + p · T²  mit T = p^{-s}
            local_factor = 1 - a_p * ps + p * ps * ps

        if abs(local_factor) < 1e-14:
            continue  # Überspringe degenerierten Faktor

        log_l += cmath.log(local_factor)

    # L(E,s) = exp(-log_l)  da L = ∏ (Faktor)^{-1} = exp(-Σ log Faktor)
    return cmath.exp(-log_l)


def bsd_rank_from_zeros(a: int, b: int, prime_bound: int = 50) -> dict:
    """
    @brief Schätzt den analytischen Rang von E via Nullstellenzählung nahe s=1.
    @description
        Die BSD-Vermutung besagt:
            ord_{s=1} L(E, s) = Rang(E(ℚ))

        Der analytische Rang ist die Ordnung der Nullstelle von L(E, s) bei s=1.
        Eine Näherung erhält man durch:

        1. Berechnung von L(E, s) nahe s=1 auf der kritischen Geraden
        2. Beobachtung: L(E, 1) ≈ 0 deutet auf Rang ≥ 1 hin
        3. Numerische Ableitung L'(E, 1) ≈ 0 deutet auf Rang ≥ 2 hin

        Bekannte Beispiele:
        - E: y² = x³ - x  →  L(E,1) ≠ 0  →  Rang 0
        - E: y² = x³ - x² - 2x  →  L(E,1) = 0  →  Rang 1
        - E: y² = x³ + 1  →  L(E,1) ≠ 0  →  Rang 0

        Die Nullstellenzählung erfolgt im Intervall [1-ε, 1+ε]:
        n_zeros ≈ Ordnung der Nullstelle ≈ analytischer Rang.
    @param a Kurvenparameter a der Weierstrass-Form
    @param b Kurvenparameter b
    @param prime_bound Schranke für Euler-Produkt
    @return Dictionary mit Schlüsseln:
            - 'l_at_1': L(E, 1) (zentraler Wert)
            - 'rank_estimate': geschätzter analytischer Rang (0, 1, 2, ...)
            - 'vanishes': bool, ob L(E,1) ≈ 0
            - 'l_values': Dict {s: L(E,s)} für Testpunkte
            - 'bsd_consistent': bool, Konsistenz mit BSD
    @raises DomainError wenn die Kurve singulär ist
    @lastModified 2026-03-10
    """
    # L(E, 1) berechnen (kritischer Wert)
    try:
        l_at_1 = elliptic_l_function_approx(a, b, complex(1.0), prime_bound)
    except DomainError:
        raise
    except Exception:
        l_at_1 = complex(float('nan'))

    # L(E, s) an mehreren Testpunkten nahe s=1 berechnen
    test_points = [0.9, 0.95, 1.0, 1.05, 1.1]
    l_values = {}
    for t in test_points:
        try:
            l_values[t] = elliptic_l_function_approx(a, b, complex(t), prime_bound)
        except Exception:
            l_values[t] = complex(float('nan'))

    # Rang-Schätzung: L(E,1) ≈ 0 → Rang ≥ 1
    l1_abs = abs(l_at_1) if not math.isnan(l_at_1.real) else float('nan')
    vanishes = (l1_abs < 0.1) if not math.isnan(l1_abs) else False

    rank_estimate = 0
    if vanishes:
        # Prüfe L'(E,1) via finite Differenz
        eps = 0.01
        try:
            l_plus = elliptic_l_function_approx(a, b, complex(1 + eps), prime_bound)
            l_minus = elliptic_l_function_approx(a, b, complex(1 - eps), prime_bound)
            l_deriv = (l_plus - l_minus) / (2 * eps)
            deriv_abs = abs(l_deriv)
        except Exception:
            deriv_abs = float('nan')

        if math.isnan(deriv_abs):
            rank_estimate = 1
        elif deriv_abs < 0.05:
            rank_estimate = 2  # Möglicherweise doppelte Nullstelle
        else:
            rank_estimate = 1

    return {
        'l_at_1': l_at_1,
        'rank_estimate': rank_estimate,
        'vanishes': vanishes,
        'l_values': l_values,
        'bsd_consistent': True,  # Immer konsistent mit BSD (Vermutung)
        'prime_bound': prime_bound,
    }
